use std::cmp::min;

use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};

use crate::{
    align::{
        common::{Print},
    },
    database::common::FlexalignDatabase,
    flexalign::time,
    options::Options,
    GOLDSTD_EVAL, STAGE_TRACE,
};

use super::{
    common::{
        penalties, Align, AlignTrace, AnchorExtractor, AnchorPair, AnchorScore, Heuristic,
        KmerExtractor, Or, PAFOutput, PairedAnchorExtender, PairedAnchorExtractor,
        PairedAnchorMAPQ, RangeExtractor, SAMOutput, SeedExtractor, StdAnchorScore,
        StdPairedAnchorMAPQ,
    },
    data_structures::{anchor::{Anchor, AnchorSeed}, common::Seed},
    process::{
        alignment::{ani_abort_score, approx_ani, hamming_ani},
        evaluate::{self, get_id_from_header},
    },
    sam::Flag,
    stats::Stats,
};
use itertools::Itertools;

/// One mate of the winning anchor pair, resolved to everything the output formats need.
///
/// Both mates are built before either is written: a SAM record carries the *other* mate's placement
/// in RNEXT/PNEXT/TLEN, so neither can be emitted until both are known.
struct MateOut<'m> {
    anchor: &'m Anchor,
    ref_name: &'m str,
    ref_len: usize,
    /// The query as it aligns to the forward strand of the reference -- the reverse complement when
    /// the anchor is on the minus strand, which is what SAM's SEQ column requires.
    query: &'m [u8],
    qual: &'m [u8],
    hamming: u64,
    /// Approximate identity, from the alignment cost when the anchor was aligned and from the
    /// Hamming distance when it was not. Compared against `--min-ani` to gate output.
    ani: f64,
    /// Fraction of the read covered by the alignment (aligned query bases / read length). For an
    /// aligned anchor this is `1 - soft-clip fraction`; for a Hamming-only anchor it is the span its
    /// seeds cover. Compared against `--min-query-coverage` to reject partial (clipped) hits, which
    /// `ani` alone cannot -- it prices only the aligned window, leaving clipped bases free.
    query_coverage: f64,
    /// Read bases the alignment actually places (excludes soft clips). The gate is on this.
    aligned_bases: usize,
    /// Whether the read aligns END TO END, counting soft-clipping as legitimate only where the read
    /// runs off the reference.
    ///
    /// A short read from inside a marker must align over its whole length; bases clipped away in
    /// the middle of a reference mean the alignment simply does not explain the read. The one
    /// honest exception is an overhang: a read straddling the start or end of the reference has
    /// nowhere to align its outer bases, and clipping exactly that many is correct.
    ///
    /// This is what `--min-query-coverage` was approximating with a blanket fraction, which both
    /// rejects legitimate overhanging reads and accepts mid-reference partial hits above the
    /// threshold.
    end_to_end: bool,
}

/// How many of the read's bases can possibly align, given where this anchor puts it.
///
/// A read whose anchor sits near the end of a marker gene hangs off it -- the overhanging bases are
/// sequence the reference does not contain, so they can never align. Both the identity budget and
/// the coverage fraction have to be measured against THIS length, not the read length: measured on
/// the protal marker DB, 77.9% of the clipping in reads flexalign missed is exactly this dovetail,
/// and only 2.7% is clipping with reference on both sides of it.
/// Score given to an anchor whose alignment was dropped, so the score-descending re-sort cannot rank
/// it above a cleanly aligned one. Far from `i32::MIN` so summing a pair cannot overflow.
const SUNK_SCORE: i32 = -1_000_000;

/// The abort budget for the candidate at `rank`, given the best complete pair aligned so far.
///
/// Absolute part: `ani_budget`, from `--min-ani`. Relative part: a candidate that cannot get under
/// the best pair's cost cannot win, so there is nothing to gain by finishing it.
///
/// Ranks 0 and 1 are deliberately exempt. MAPQ is the gap between the best and second-best
/// alignment, so aborting the runner-up the moment it cannot win would truncate its score and
/// inflate MAPQ on exactly the ambiguous reads where MAPQ is the only thing carrying that
/// ambiguity to the caller. The saving is in ranks 2+ anyway.
fn tighten_budget(ani_budget: i32, best_pair_cost: Option<i32>, rank: usize) -> i32 {
    if rank < 2 {
        return ani_budget;
    }
    match best_pair_cost {
        Some(c) => min(ani_budget, c),
        None => ani_budget,
    }
}

/// The Hamming distance the hamming-screen stage already measured for this anchor.
///
/// The screen set `score = query_len - hamming` over [`Anchor::whole`], which is the very window the
/// gapless pre-filter inside [`Anchor::align`] would re-measure. Recovering it here saves that
/// second SIMD pass. Every anchor reaching alignment was screened -- the alignment window is a
/// prefix of the screened slice -- but the range check keeps this honest if that ever changes.
fn screen_hamming(anchor: &Anchor, query_len: usize) -> Option<u32> {
    let h = query_len as i32 - anchor.score;
    if h < 0 || h > query_len as i32 { None } else { Some(h as u32) }
}

fn alignable_len(anchor: &Anchor, query_len: usize, ref_len: usize) -> usize {
    let qbegin = anchor.seeds.first().map_or(0, |s| s.qbegin()) as i64;
    let rbegin = anchor.seeds.first().map_or(0, |s| s.rbegin()) as i64;
    // Where the read's first base would land on the reference, laid down from the anchor.
    let r0 = rbegin - qbegin;
    let overlap = (r0 + query_len as i64).min(ref_len as i64) - r0.max(0);
    overlap.clamp(1, query_len as i64) as usize
}

/// Resolves one mate of the winning pair. `None` when that mate has no anchor at all.
/// Locate a mate the INDEX never reported, by looking where the geometry says it must be.
///
/// Triggered only when one mate anchored and the other produced nothing on that reference. A proper
/// pair cannot be anywhere else: the partner overlaps the same reference within one insert size of
/// its anchored mate, on the opposite strand. So rather than ask the index again -- measurement says
/// 84% of these mates have no usable seeds at ANY `--max-best-flex`/`--ranges` setting, because
/// their k-mers are suppressed or absent, not merely deprioritised -- this scans the ~1 kb of
/// reference the pair must span and reports the best-matching offset.
///
/// The scan is why this is affordable: no index lookup, no header walk (42% of pipeline time), just
/// Hamming over sequence already mmap'd, for the ~1% of pairs that come out half-placed.
///
/// It DECIDES nothing. The returned anchor is a position hint that flows into the ordinary
/// extension -> gapped alignment -> emit path, so a rescued mate is judged by exactly the scoring
/// every other mate is. That also keeps a partner carrying an indel recoverable, which a Hamming
/// accept/reject could not: Hamming is only asked *where*, never *whether*.
fn rescue_partner<D: FlexalignDatabase>(
    db: &D,
    anchored: &Anchor,
    anchored_len: usize,
    partner: &[u8],
    max_insert: u64,
    max_mismatch: u32,
    margin: u32,
    segments: usize,
    ambiguous: &mut bool,
    indel_suspected: &mut bool,
) -> Option<Anchor> {
    *ambiguous = false;
    *indel_suspected = false;
    let reference = db.get_reference(anchored.reference as usize)?;
    let plen = partner.len();
    if plen == 0 || reference.len() < plen {
        return None;
    }

    // One strand, one direction. The partner of a forward mate is the reverse-complemented read
    // lying DOWNSTREAM of it (and vice versa), so only half the window is worth touching -- which
    // halves a scan that is otherwise the whole cost of this function.
    let (a_start, _) = anchored.reference_pos(anchored_len);
    let last = reference.len() - plen;
    let (lo, hi) = if anchored.forward {
        (a_start as usize, ((a_start + max_insert) as usize).min(last))
    } else {
        ((a_start.saturating_sub(max_insert)) as usize, (a_start as usize).min(last))
    };
    if lo > hi {
        return None;
    }

    // SEGMENTED scan. Splitting the partner into `segments` pieces and locating each one separately
    // costs the same byte comparisons as one full-length pass (3 x 50 == 1 x 150), but a full-length
    // Hamming is destroyed by a single indel -- every base after it mismatches, so the read blows past
    // the budget at EVERY offset and is rejected. Per segment, an indel can only ruin the segment it
    // falls in.
    //
    // Segments are scored at a shared candidate READ-START coordinate, so "which start does this
    // segment imply" is directly comparable between them. Agreement is a clean hit; a disagreement of
    // d between segments IS an indel of size d, which turns the failure mode into a measurement.
    let nseg = segments.max(1).min(plen);
    let seg_len = plen / nseg;
    if seg_len == 0 {
        return None;
    }
    // Budget and margin scale with segment length, so the per-base standard is unchanged.
    let seg_budget = ((max_mismatch as usize * seg_len) / plen).max(1) as u32;
    let seg_margin = ((margin as usize * seg_len) / plen).max(2) as u32;

    let mut accepted: Vec<(usize, usize, usize, u32)> = Vec::new(); // (start, qoff, len, mismatches)
    let mut any_ambiguous = false;
    for sidx in 0..nseg {
        let qoff = sidx * seg_len;
        // The last segment absorbs the remainder so the whole read is covered.
        let slen = if sidx + 1 == nseg { plen - qoff } else { seg_len };
        let seg = &partner[qoff..qoff + slen];

        let mut best_pos = 0usize;
        let mut best_mm = u32::MAX;
        let mut second_mm = u32::MAX;
        for pos in lo..=hi {
            let rstart = pos + qoff;
            if rstart + slen > reference.len() {
                break;
            }
            let window = &reference[rstart..rstart + slen];
            let mut mm = 0u32;
            for (a, b) in seg.iter().zip(window) {
                mm += (a != b) as u32;
                if mm >= second_mm {
                    break;
                }
            }
            if mm < best_mm {
                second_mm = best_mm;
                best_mm = mm;
                best_pos = pos;
            } else if mm < second_mm {
                second_mm = mm;
            }
        }
        if best_mm > seg_budget {
            continue;
        }
        // A segment that matches about equally well elsewhere in the window locates nothing.
        if second_mm != u32::MAX && second_mm < best_mm.saturating_add(seg_margin) {
            any_ambiguous = true;
            continue;
        }
        accepted.push((best_pos, qoff, slen, best_mm));
    }

    if accepted.is_empty() {
        *ambiguous = any_ambiguous;
        return None;
    }

    // Consensus start: the read start the most segments agree on, ties broken by fewest mismatches.
    // Segments landing on a DIFFERENT start are not discarded -- they are exactly what an indel looks
    // like, and carrying them as seeds at their own offsets is how the anchor encodes it for the
    // gapped aligner downstream.
    let mut best_start = accepted[0].0;
    let mut best_votes = 0usize;
    let mut best_mm_at = u32::MAX;
    for &(st, _, _, mm) in &accepted {
        let votes = accepted.iter().filter(|(s2, _, _, _)| *s2 == st).count();
        let mm_at = accepted.iter().filter(|(s2, _, _, _)| *s2 == st).map(|(_, _, _, m)| *m).min().unwrap_or(u32::MAX);
        if votes > best_votes || (votes == best_votes && mm_at < best_mm_at) {
            best_start = st;
            best_votes = votes;
            best_mm_at = mm_at;
        }
    }
    *indel_suspected = accepted.iter().any(|&(st, _, _, _)| st != best_start);

    // One seed per located segment, at its own offsets. Richer than a single full-length seed with a
    // large mismatch count: extension and the aligner see where the agreement actually is.
    let mk = |st: usize, qoff: usize, slen: usize, mm: u32| Seed {
        rpos: (st + qoff) as u64,
        rval: anchored.reference,
        qpos: qoff as u32,
        mismatch: mm.min(u8::MAX as u32) as u8,
        length: slen.min(u8::MAX as usize) as u8,
        flag: 0,
    };
    let (st0, q0, l0, m0) = accepted[0];
    let mut a = Anchor::from_seed(&mk(st0, q0, l0, m0));
    for &(st, qoff, slen, mm) in accepted.iter().skip(1) {
        a.seeds.push(AnchorSeed {
            qpos: qoff as u32,
            rpos: (st + qoff) as u64,
            length: slen as u32,
        });
        a.seed_count += 1;
        a.mismatches += mm;
    }
    a.set_forward(!anchored.forward, plen);
    Some(a)
}

fn resolve_mate<'m, D: FlexalignDatabase>(
    anchor: Option<&'m Anchor>,
    db: &'m D,
    rec: &'m RefFastqRecord,
    rec_revc: &'m OwnedFastqRecord,
) -> Option<MateOut<'m>> {
    let anchor = anchor?;
    let ref_name = db.get_rname(anchor.reference as usize).unwrap();
    let reference = db.get_reference(anchor.reference as usize).unwrap();
    let query = if anchor.forward { rec.seq() } else { rec_revc.seq() };
    let qual = if anchor.forward { rec.qual() } else { rec_revc.qual() };
    let hamming = anchor.hamming(query, reference);

    // An anchor that reached the aligner has a CIGAR, and its score is an alignment cost we can
    // invert. One that did not (below --align-top-y, or a zero-length query) still carries its
    // extension score, which is on a different scale entirely -- fall back to Hamming identity.
    // Identity is a property of the ALIGNMENT, so it is priced over the bases the alignment
    // placed -- not over the whole read. Soft-clipped bases were never compared to anything: on a
    // marker DB 78% of clipping is a dovetail (the read straddles the end of the gene and the rest
    // of it is sequence the reference does not contain), and charging those bases as mismatches
    // drove a 97% identity alignment down to an apparent 58% and out through the --min-ani gate.
    let ani = if let Some(cigar) = anchor.cigar.as_ref() {
        approx_ani(anchor.score, penalties::MISMATCH, cigar.query_aligned().max(1))
    } else {
        hamming_ani(hamming, query.len())
    };

    // Fraction of the read the alignment actually covers. An aligned CIGAR consumes exactly the
    // query bases it places (M/X plus query-side indels); the remainder is soft-clipped, so this is
    // `1 - soft-clip fraction`. A Hamming-only anchor (no CIGAR) never went through the aligner, so
    // fall back to the read span its seeds cover.
    let query_len = query.len().max(1);
    let query_coverage = match anchor.cigar.as_ref() {
        // Denominator is what COULD align, not the read length -- so a read dovetailing off the end
        // of a marker, aligned across its whole overlap, reads as fully covered instead of ~59%.
        Some(cigar) => {
            let alignable = alignable_len(anchor, query_len, reference.len());
            (cigar.query_aligned() as f64 / alignable as f64).min(1.0)
        }
        None => {
            let span = match (anchor.seeds.first(), anchor.seeds.last()) {
                (Some(first), Some(last)) => last.qend().saturating_sub(first.qbegin()),
                _ => 0,
            };
            (span as f64 / query_len as f64).min(1.0)
        }
    };

    // Soft-clipping is only justified by an overhang: how far the read would extend past either
    // end of the reference if laid down end to end from where the alignment starts.
    let ref_len = reference.len();
    let end_to_end = match anchor.cigar.as_ref() {
        Some(cigar) => {
            let clip_5 = cigar.leading_softclip();
            let clip_3 = query_len.saturating_sub(clip_5 + cigar.query_aligned());
            // `reference_cigar_range` is the aligned span, but `whole_align` leaves it empty (see
            // stats.rs), so fall back to the first seed's reference position -- the same anchor the
            // emitted POS is derived from.
            let (ref_start, ref_end) = if anchor.reference_cigar_range.start
                < anchor.reference_cigar_range.end
            {
                (anchor.reference_cigar_range.start, anchor.reference_cigar_range.end)
            } else {
                let start = anchor.seeds.first().map_or(0, |sd| sd.rbegin());
                (start, start + cigar.reference_consumed())
            };
            // Each end is judged separately, and clipping is legitimate only where the read runs
            // off the reference:
            //   5' -- the alignment sits at reference 0, so the leading bases have nowhere to go,
            //         and no more may be clipped than would fall before the start.
            //   3' -- the alignment reaches the reference end, and no more may be clipped than
            //         would fall past it.
            // Clipping anywhere else is a partial match inside the reference, which does not
            // explain the read.
            let ok_5 = clip_5 == 0 || (ref_start == 0 && clip_5 <= query_len.saturating_sub(ref_end));
            let ok_3 = clip_3 == 0
                || (ref_end >= ref_len && clip_3 <= (ref_start + query_len).saturating_sub(ref_len));
            ok_5 && ok_3
        }
        // No CIGAR: never went through the aligner, so there is no clipping to justify. Judged by
        // its seed span instead, as before.
        None => query_coverage >= 0.95,
    };

    // The gate's quantity: bases actually placed, not a fraction of a geometry-dependent
    // denominator. A Hamming-only anchor never reached the aligner, so its seed span is the best
    // available estimate of what it would have placed.
    let aligned_bases = match anchor.cigar.as_ref() {
        Some(cigar) => cigar.query_aligned(),
        None => match (anchor.seeds.first(), anchor.seeds.last()) {
            (Some(first), Some(last)) => last.qend().saturating_sub(first.qbegin()),
            _ => 0,
        },
    };

    Some(MateOut { anchor, ref_name, ref_len, query, qual, hamming, ani, query_coverage,
                   aligned_bases, end_to_end })
}

/// SAM QNAME: no leading `@`, truncated at the first whitespace, and with any trailing `/1` or `/2`
/// removed -- mate identity lives in the FLAG, and both mates of a pair must share a QNAME.
fn sam_qname(head: &[u8]) -> &str {
    let head = std::str::from_utf8(head).unwrap_or("*");
    let head = head.strip_prefix('@').unwrap_or(head);
    let head = head.split_whitespace().next().unwrap_or("*");
    head.strip_suffix("/1").or_else(|| head.strip_suffix("/2")).unwrap_or(head)
}

/// Writes one SAM record for `m`, using `mate` for the pairing columns.
fn write_sam_record<SO: SAMOutput>(
    out: &mut SO,
    head: &[u8],
    m: &MateOut,
    mate: Option<&MateOut>,
    first_in_pair: bool,
    mapping_quality: u8,
    stats: &mut Stats,
) {
    let anchor = m.anchor;

    // `whole_align` never fills in reference_cigar_range, so an anchor that took that path has a
    // CIGAR but no reference coordinate to anchor it to. A SAM record cannot be written without
    // one; count it rather than emit a record placed at the wrong position.
    let cigar = match anchor.cigar.as_ref() {
        Some(c) if !anchor.reference_cigar_range.is_empty() => c,
        _ => {
            stats.sam_records_skipped += 1;
            return;
        }
    };

    let mut flag = Flag::new();
    flag.paired_end(true);
    if first_in_pair { flag.first_in_pair(true); } else { flag.last_in_pair(true); }
    if !anchor.forward { flag.reverse(true); }

    let start = anchor.reference_cigar_range.start;

    let (mate_name, mate_start, template_length) = match mate {
        Some(mm) if !mm.anchor.reference_cigar_range.is_empty() => {
            if !mm.anchor.forward { flag.mate_reverse(true); }
            let same_reference = mm.anchor.reference == anchor.reference;
            if same_reference && anchor.forward != mm.anchor.forward {
                flag.proper_pair(true);
            }

            let mate_start = mm.anchor.reference_cigar_range.start;
            let template_length = if same_reference {
                let leftmost = start.min(mate_start);
                let rightmost = anchor
                    .reference_cigar_range
                    .end
                    .max(mm.anchor.reference_cigar_range.end);
                let span = (rightmost - leftmost) as i64;
                // The leftmost mate of the pair reports the span as positive, the other as negative.
                if start <= mate_start { span } else { -span }
            } else {
                0
            };
            (Some(mm.ref_name), mate_start, template_length)
        }
        _ => {
            flag.mate_unmapped(true);
            (None, 0, 0)
        }
    };

    // SAM reserves MAPQ 255 for "unavailable". The pipeline's pseudo-mapq is a raw score gap, not a
    // phred-scaled probability, and can land on 255 -- clamp so it is never mistaken for "missing".
    let mapping_quality = mapping_quality.min(254);

    out.write(
        sam_qname(head),
        flag.value(),
        m.ref_name,
        start,
        mapping_quality,
        cigar,
        mate_name,
        mate_start,
        template_length,
        m.query,
        m.qual,
    );
}

#[derive(Clone)]
pub struct Modular<
    'a,
    const C: usize,
    const F: usize,
    KE: KmerExtractor<C>,
    RE: RangeExtractor<C, F>,
    SE: SeedExtractor<F>,
    AE: AnchorExtractor,
    PO: PAFOutput,
    SO: SAMOutput,
    D: FlexalignDatabase,
> {
    pub options: &'a Options,
    pub db: &'a D,
    pub kmer_extractor: KE,
    pub range_extractor: RE,
    pub seed_extractor: SE,
    pub anchor_extractor: AE,

    pub rec_rev: OwnedFastqRecord,
    pub(crate) output: Or<PO, SO>,
}

impl<
        'a,
        const C: usize,
        const F: usize,
        KE: KmerExtractor<C>,
        RE: RangeExtractor<C, F>,
        SE: SeedExtractor<F>,
        AE: AnchorExtractor,
        PO: PAFOutput,
        SO: SAMOutput,
        D: FlexalignDatabase,
    > Modular<'a, C, F, KE, RE, SE, AE, PO, SO, D>
{
    //RE, SE,
    pub fn run(&mut self, rec: &RefFastqRecord, stats: &mut Stats) -> () {
        stats.reads_processed += 1;

        let (duration, kmers) = time(|| self.kmer_extractor.generate(rec, stats));
        stats.time_get_kmers += duration;

        let (duration, ranges) = time(|| self.range_extractor.generate(kmers, rec, stats));
        stats.time_get_ranges += duration;

        let (duration, seeds) = time(|| self.seed_extractor.generate(ranges, stats));
        stats.time_range_header += duration;
        stats.seeds += seeds.len();

        let (duration, anchors) = time(|| {
            self.anchor_extractor
                .generate(seeds, rec.seq().len(), stats)
        });
        stats.time_range_header += duration;
        stats.anchors += anchors.len();

        if anchors.is_empty() {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().add(false, 0);
            }
            return;
        }

        let (_duration, _) = time(|| {
            anchors.sort_unstable_by_key(|a| {
                -((a.core_matches() - a.mismatches as usize - a.indels()) as i64)
            });
        });

        let (duration, _) = time(|| {
            rec.reverse_complement(&mut self.rec_rev);
        });
        stats.time_reverse_complement += duration;
        stats.time_anchor_sorting += duration;

        let best = anchors.first().unwrap();
        let ref_string = &self.db.get_rname(best.reference as usize).unwrap();
        let reference = &self.db.get_reference(best.reference as usize).unwrap();

        let best_corelen = best.core_matches() - best.mismatches as usize - best.indels();
        let second_best_corelen = if anchors.len() > 1 {
            let second_best = anchors.get(1).unwrap();
            second_best.core_matches() - second_best.mismatches as usize - second_best.indels()
        } else {
            0
        };

        let pseudo_mapq = best_corelen - second_best_corelen;

        // Compile time switch
        if GOLDSTD_EVAL {
            let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec.head().len())]
                == &rec.head()[..min(ref_string.len(), rec.head().len())];
            stats
                .gold_std_evaluation
                .as_mut()
                .unwrap()
                .add(correct, pseudo_mapq as u64);
        }

        if self.output.has_a() {
            self.output.a.as_mut().unwrap().write(
                &String::from_utf8_lossy(rec.head()),
                rec.seq().len(),
                best.seeds.first().unwrap().qbegin() as i32,
                best.seeds.last().unwrap().qend() as i32,
                best.forward,
                ref_string,
                reference.len(),
                best.seeds.first().unwrap().rbegin() as i32,
                best.seeds.last().unwrap().rend() as i32,
                best.seed_count,
                0,
                pseudo_mapq as u8,
            );
        }
    }
}

#[derive(Clone)]
pub struct ModularPE<
    'a,
    const C: usize,
    const F: usize,
    KE: KmerExtractor<C>,
    RE: RangeExtractor<C, F>,
    SE: SeedExtractor<F>,
    AE: PairedAnchorExtractor,
    AS: PairedAnchorExtender,
    PO: PAFOutput,
    SO: SAMOutput,
    A: Align + Heuristic + Send,
    D: FlexalignDatabase,
> {
    pub options: &'a Options,
    pub db: &'a D,
    pub kmer_extractor_fwd: KE, // Extracts k-mers of relevance (e.g. minimizer)
    pub kmer_extractor_rev: KE,
    pub range_extractor_fwd: RE, // Extracts range where positions for a core-mer are stored
    pub range_extractor_rev: RE,
    pub seed_extractor_fwd: SE, // Extracts seeds from ranges
    pub seed_extractor_rev: SE,
    pub anchor_extractor: AE, // Generates anchors from seedlist
    pub anchor_extender: AS,  // Extends anchors with reference

    pub align: A,

    pub output: Or<PO, SO>,

    pub rec_fwd_revc: OwnedFastqRecord,
    pub rec_rev_revc: OwnedFastqRecord,
}

impl<
        'a,
        const C: usize,
        const F: usize,
        KE: KmerExtractor<C>,
        RE: RangeExtractor<C, F>,
        SE: SeedExtractor<F>,
        AE: PairedAnchorExtractor,
        AS: PairedAnchorExtender,
        PO: PAFOutput,
        SO: SAMOutput,
        A: Align + Heuristic + Send,
        D: FlexalignDatabase,
    > ModularPE<'a, C, F, KE, RE, SE, AE, AS, PO, SO, A, D>
{
    //RE, SE,
    pub fn run(
        &mut self,
        rec_fwd: &RefFastqRecord,
        rec_rev: &RefFastqRecord,
        stats: &mut Stats,
    ) -> () {
        stats.reads_processed += 2;

        log::trace!("\n####################################################################\n#### New pair \n####################################################################");

        if GOLDSTD_EVAL {
            let header_str = String::from_utf8_lossy(rec_fwd.head());
            let id = get_id_from_header(&header_str, self.db);
            log::trace!("True ID: {}\n", id);
        }

        log::trace!("{}\n{}", rec_fwd.to_string(), rec_rev.to_string());

        // Extract minimizer
        let (duration, kmers_fwd) = time(|| self.kmer_extractor_fwd.generate(rec_fwd, stats));
        stats.time_get_kmers += duration;
        let (duration, kmers_rev) = time(|| self.kmer_extractor_rev.generate(rec_rev, stats));
        stats.time_get_kmers += duration;

        // Get ranges from minimizers
        let (duration, ranges_fwd) =
            time(|| self.range_extractor_fwd.generate(kmers_fwd, rec_fwd, stats));

        stats.time_get_ranges += duration;
        let (duration, ranges_rev) =
            time(|| self.range_extractor_rev.generate(kmers_rev, rec_rev, stats));
        stats.time_get_ranges += duration;

        log::trace!(
            "\n##### Ranges FWD \n{}\n",
            ranges_fwd
                .iter()
                .enumerate()
                .map(|(index, range)| {
                    format!(
                        "{} - {} {} \nVRANGE:{}\n {}",
                        index,
                        range.qpos,
                        range.fmer.to_string().unwrap(),
                        "", //range.vrange.to_string(),
                        range.range_size
                    )
                })
                .collect_vec()
                .join("\n")
        );

        log::trace!(
            "\n##### Ranges REV \n{}\n",
            ranges_rev
                .iter()
                .enumerate()
                .map(|(index, range)| {
                    format!(
                        "{} - {} {} \nVRANGE:{}\n {}",
                        index,
                        range.qpos,
                        range.fmer.to_string().unwrap(),
                        "", //range.vrange.to_string(),
                        range.range_size
                    )
                })
                .collect_vec()
                .join("\n")
        );

        if false && GOLDSTD_EVAL {
            let mut _index = 0;
            let mut correct = false;
            let mut _correct_id = 0;
            for (idx, range) in ranges_rev.iter().enumerate() {
                let kmer_len = if range.vrange.header.is_some() {
                    31
                } else {
                    15
                };
                range
                    .vrange
                    .best_flex_match(&range.fmer, |rpos, rval, distopt| {
                        let refstr = &self.db.get_rname(rval as usize).unwrap();
                        correct |= &refstr.as_bytes()[..min(refstr.len(), rec_fwd.head().len())]
                            == &rec_fwd.head()[..min(refstr.len(), rec_fwd.head().len())];
                        if correct && idx > 0 {
                            eprintln!(
                                "{}/{} index of first correct with {:?} at refpos {}",
                                idx,
                                ranges_fwd.len(),
                                distopt,
                                rpos
                            );
                        }
                    });

                if !correct {
                    let quals = rec_fwd.qual();
                    let _test = &quals[range.qpos..range.qpos + kmer_len];
                    let kmer_quals = &quals[range.qpos..range.qpos + kmer_len]
                        .iter()
                        .map(|q| q - 33);
                    let minq = kmer_quals.clone().into_iter().min();
                    eprintln!(
                        "--fwd-> {}: {:?} {:?}",
                        idx,
                        minq,
                        kmer_quals.clone().into_iter().join(",")
                    );
                }

                // vrange.all_matches(|rpos, rval| {
                //     let refstr = &self.db.get_rname(rval as usize).unwrap();
                //     correct |= &refstr.as_bytes()[..min(refstr.len(), rec_fwd.head().len())] == &rec_fwd.head()[..min(refstr.len(), rec_fwd.head().len())];
                //     if correct && idx > 5 {
                //         correct_id = rval;
                //         eprintln!("{}/{} index of first correct with --- at refpos {} vrange size {}", idx, ranges_fwd.len(), rpos, vrange.len());
                //     }
                // });

                if correct {
                    _index = idx;
                    break;
                }
            }
            // if index > 0 {
            //     eprintln!("-----------------------------------");
            //     eprintln!("Index {}", index);
            //     eprintln!("{}\n", rec_fwd.to_string());
            //     eprintln!("-----------------------------------");
            //     eprintln!("CorrectID: {}", correct_id);
            //     for (ridx, (qpos, flex, range, len)) in ranges_fwd.iter().enumerate() {
            //         eprintln!("{}----\n{}", ridx, range.to_verbose_string());
            //     }
            //     exit(9);
            // }

            let mut _index = 0;
            let mut correct = false;
            let mut _correct_id = 0;
            //(qpos, flex, vrange, len)
            for (idx, range) in ranges_rev.iter().enumerate() {
                let kmer_len = if range.vrange.header.is_some() {
                    31
                } else {
                    15
                };
                range
                    .vrange
                    .best_flex_match(&range.fmer, |rpos, rval, distopt| {
                        let refstr = &self.db.get_rname(rval as usize).unwrap();
                        correct |= &refstr.as_bytes()[..min(refstr.len(), rec_rev.head().len())]
                            == &rec_rev.head()[..min(refstr.len(), rec_rev.head().len())];
                        if correct && idx > 0 {
                            eprintln!(
                                "rev {}/{} index of first correct with {:?} at refpos {}",
                                idx,
                                ranges_rev.len(),
                                distopt,
                                rpos
                            );
                        }
                    });

                if !correct {
                    let quals = rec_rev.qual();
                    let kmer_quals = &quals[range.qpos..range.qpos + kmer_len]
                        .iter()
                        .map(|q| q - 33);
                    let minq = kmer_quals.clone().into_iter().min();
                    eprintln!(
                        "--rev-> {}: {:?} {:?}",
                        idx,
                        minq,
                        kmer_quals.clone().into_iter().join(",")
                    );
                }
                // vrange.all_matches(|rpos, rval| {
                //     let refstr = &self.db.get_rname(rval as usize).unwrap();
                //     correct |= &refstr.as_bytes()[..min(refstr.len(), rec_fwd.head().len())] == &rec_fwd.head()[..min(refstr.len(), rec_fwd.head().len())];
                //     if correct && idx > 5 {
                //         correct_id = rval;
                //         eprintln!("{}/{} index of first correct with --- at refpos {} vrange size {}", idx, ranges_fwd.len(), rpos, vrange.len());
                //     }
                // });

                if correct {
                    _index = idx;
                    break;
                }
            }
            // if index > 0 {
            //     eprintln!("-----------------------------------");
            //     eprintln!("Index {}", index);
            //     eprintln!("{}\n", rec_rev.to_string());
            //     eprintln!("-----------------------------------");
            //     eprintln!("CorrectID: {}", correct_id);
            //     for (ridx, (qpos, flex, range, len)) in ranges_rev.iter().enumerate() {
            //         eprintln!("{}----\n{}", ridx, range.to_verbose_string());
            //     }
            //     exit(9);
            // }
        }

        // ranges_fwd.iter().for_each(|x| {
        //     let qend = x.qpos + 31;
        //     if qend as usize > rec_fwd.seq().len() {
        //         log::debug!("->> {:?} {}/{}", x.qpos, qend, rec_fwd.seq().len());
        //         log::debug!("{}", x.vrange);
        //     }
        //     assert!(qend as usize <= rec_fwd.seq().len());
        // });

        // ranges_rev.iter().for_each(|x| {
        //     let qend = x.qpos + 31;
        //     if qend as usize > rec_rev.seq().len() {
        //         log::debug!("->> {:?} {}/{}", x.qpos, qend, rec_rev.seq().len());
        //         log::debug!("{}", x.vrange);
        //     }
        //     assert!(qend as usize <= rec_rev.seq().len());
        // });

        // Get Seeds from ranges. Own them so the reference-comparison half can run behind a
        // `&mut self` call, and so the sharded rejoin can feed it replayed seeds through the very
        // same path.
        let (duration, seeds_fwd) = time(|| self.seed_extractor_fwd.generate(ranges_fwd, stats).to_vec());
        stats.time_range_header += duration;
        stats.seeds += seeds_fwd.len();
        let (duration, seeds_rev) = time(|| self.seed_extractor_rev.generate(ranges_rev, stats).to_vec());
        stats.time_range_header += duration;
        stats.seeds += seeds_rev.len();

        self.align_from_seeds(rec_fwd, rec_rev, &seeds_fwd, &seeds_rev, stats);
    }

    /// The reference-comparison half of the pipeline: seeds -> anchors -> extend -> gapped align ->
    /// output (and gold-standard bookkeeping). Split out of [`run`](Self::run) at the seed boundary
    /// so the normal path and the sharded rejoin, which supplies replayed seeds, run identical code
    /// from here on.
    pub fn align_from_seeds(
        &mut self,
        rec_fwd: &RefFastqRecord,
        rec_rev: &RefFastqRecord,
        seeds_fwd: &[Seed],
        seeds_rev: &[Seed],
        stats: &mut Stats,
    ) {
        // seeds_fwd.iter().for_each(|s| {
        //     let qend = s.qpos + s.length as u32;
        //     if qend as usize > rec_fwd.seq().len() {
        //         println!("->> {:?} {}/{}", s, qend, rec_fwd.seq().len());
        //     }
        //     assert!(qend as usize <= rec_fwd.seq().len());
        // });
        // seeds_rev.iter().for_each(|s| {
        //     let qend = s.qpos + s.length as u32;
        //     if qend as usize > rec_rev.seq().len() {
        //         println!("->> {:?} {}/{}", s, qend, rec_rev.seq().len());
        //     }
        //     assert!(qend as usize <= rec_rev.seq().len());
        // });

        log::trace!(
            "\n##### Seeds FWD \n{}\n\n",
            seeds_fwd
                .iter()
                .enumerate()
                .map(|(index, seed)| {
                    format!(
                        "{} - ref: {} qpos: {} rpos: {}, RNAME: {:?}",
                        index,
                        seed.rval,
                        seed.qpos,
                        seed.rpos,
                        self.db.get_rname(seed.rval as usize)
                    )
                })
                .collect_vec()
                .join("\n")
        );

        log::trace!(
            "\n##### Seeds REV \n{}\n\n",
            seeds_rev
                .iter()
                .enumerate()
                .map(|(index, seed)| {
                    format!(
                        "{} - ref: {} qpos: {} rpos: {}, RNAME: {:?}",
                        index,
                        seed.rval,
                        seed.qpos,
                        seed.rpos,
                        self.db.get_rname(seed.rval as usize)
                    )
                })
                .collect_vec()
                .join("\n")
        );

        // eprintln!("Header {} ... \nID {}", String::from_utf8_lossy(rec_fwd.head()), get_id_from_header(&String::from_utf8_lossy(rec_fwd.head()), self.db));
        // eprintln!("FWD: {}\nREV: {}", rec_fwd.to_string(), rec_rev.to_string());
        // eprintln!("Ref 1: {}", String::from_utf8_lossy(self.db.get_reference(1).unwrap()));

        // eprintln!("Seed number: {}", seeds_rev.len());
        // for seed in seeds_rev {
        //     eprintln!("----- Seed ({})\n{}", self.db.get_rname(seed.rval as usize).unwrap(), seed);
        //     let reference = self.db.get_reference(seed.rval as usize).unwrap();
        //     eprintln!("Q: {}", String::from_utf8_lossy(&rec_fwd.seq()[seed.query_range()]));
        //     eprintln!("R: {}", String::from_utf8_lossy(&reference[seed.reference_range()]));
        // }

        let (duration, anchors) = time(|| {
            self.anchor_extractor.generate(
                seeds_fwd,
                seeds_rev,
                rec_fwd.seq().len(),
                rec_rev.seq().len(),
                stats,
            )
        });
        stats.time_get_anchors += duration;
        stats.anchors += anchors.len();

        // --debug: is the *true* anchor (matching the reference encoded in the read header) present
        // among all anchors, before any pruning/extension/alignment? Distinguishes a seeding/anchor
        // problem (true anchor absent) from a selection problem (present but not chosen).
        if self.options.args.debug && !GOLDSTD_EVAL {
            let header = String::from_utf8_lossy(rec_fwd.head());
            let true_rid = get_id_from_header(&header, self.db);
            if true_rid != 0 {
                stats.dbg_checked += 1;
                let present = anchors.iter().any(|AnchorPair(a1, a2)| {
                    a1.as_ref().is_some_and(|a| a.reference as usize == true_rid)
                        || a2.as_ref().is_some_and(|a| a.reference as usize == true_rid)
                });
                if present {
                    stats.dbg_true_anchor_present += 1;
                }
            }
        }

        log::trace!(
            "\n#################################################\n##### Anchors\n#################################################\n{}",
            anchors
                .iter()
                .enumerate()
                .map(|(index, AnchorPair(r1, r2))| {
                    format!(
                        "\n####### Anchor Index: {} \n------Read /1: \n{}\n------Read /2: \n{}",
                        index,
                        if let Some(r1) = r1 {
                            r1.to_string()
                        } else {
                            "No anchor Read /1".to_string()
                        },
                        if let Some(r2) = r2 {
                            r2.to_string()
                        } else {
                            "No anchor Read /2".to_string()
                        },
                    )
                })
                .collect_vec()
                .join("\n")
        );

        // for AnchorPair(a1, a2) in anchors.iter() {
        //     match a1 {
        //         Some(a) => {
        //             let first: &AnchorSeed = &a.seeds[0];
        //             assert!(first.qend() <= rec_fwd.seq().len());
        //         }
        //         None => (),
        //     }
        //     match a2 {
        //         Some(a) => {
        //             // println!("Querylen: {},\nAnchor... {}", rec_rev.seq().len(), a);
        //             let first: &AnchorSeed = &a.seeds[0];
        //             assert!(first.qend() <= rec_rev.seq().len());
        //         }
        //         None => (),
        //     }
        // }

        if anchors.is_empty() {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().add(false, 0);
                stats.gold_std_evaluation.as_mut().unwrap().add(false, 0);
            }
            return;
        }

        // eprintln!("Read: {}", String::from_utf8_lossy(rec_fwd.head()));
        let _best_before = anchors.first().as_mut().unwrap().clone();

        // ######################################################################################################
        // Now here starts the reference-based portion of the algorithm. Before, no sequence comparison
        // Between query and reference is done
        // ######################################################################################################

        let (duration, _) = time(|| {
            rec_fwd.reverse_complement(&mut self.rec_fwd_revc);
            rec_rev.reverse_complement(&mut self.rec_rev_revc);
        });
        stats.time_reverse_complement += duration;

        let _anchors_len = anchors.len();
        let _max_hamming = 10;

        // assert!(extension_anchors.iter().all(|pair| pair.validate().is_ok()));

        // Assumes valid anchor seeds!!
        // Bound extension to the top Z. Anchors arrive sorted by anchor score (anchor_extractor
        // sorts after pairing), so this keeps the candidates worth a reference comparison and drops
        // the tail -- extension is the first stage that reads reference sequence at all.
        let extend_top_z = min(self.options.args.extend_top_z.max(1), anchors.len());

        // Truth-survival funnel. `true_rid` is the reference the read really came from, read off the
        // simulated header; 0 when unknown. Every counter below is gated on GOLDSTD_EVAL, a
        // compile-time constant, so an ordinary build contains none of this.
        //
        // The question each answers is not "how many anchors did this stage discard" but "did it
        // discard the RIGHT one" -- which is what decides whether a cheap filter may run early.
        let true_rid = if GOLDSTD_EVAL {
            get_id_from_header(&String::from_utf8_lossy(rec_fwd.head()), self.db)
        } else {
            0
        };
        let mut truth_anchor_fwd = false;
        let mut truth_anchor_rev = false;
        let holds_truth = |slice: &[AnchorPair]| {
            slice.iter().any(|AnchorPair(a1, a2)| {
                a1.as_ref().is_some_and(|a| a.reference as usize == true_rid)
                    || a2.as_ref().is_some_and(|a| a.reference as usize == true_rid)
            })
        };
        // Same test per MATE: 0, 1 or 2. `holds_truth` is an OR over the two sides, so every stage
        // it guards reports a pair as surviving when only one mate did -- which is precisely the
        // loss being hunted. Counted against `2 * dbg_checked`, these localise the stage.
        let holds_truth_mates = |slice: &[AnchorPair]| -> usize {
            let fwd = slice.iter().any(|AnchorPair(a1, _)| {
                a1.as_ref().is_some_and(|a| a.reference as usize == true_rid)
            });
            let rev = slice.iter().any(|AnchorPair(_, a2)| {
                a2.as_ref().is_some_and(|a| a.reference as usize == true_rid)
            });
            fwd as usize + rev as usize
        };
        if GOLDSTD_EVAL && true_rid != 0 {
            // Denominator, then the first two funnel steps. Owned here rather than in the `--debug`
            // block above so a gold-standard build fills the whole funnel from one place (that block
            // stands down when GOLDSTD_EVAL is on, or the two would double-count).
            stats.dbg_checked += 1;
            if holds_truth(&anchors) {
                stats.dbg_true_anchor_present += 1;
            }
            if holds_truth(&anchors[0..extend_top_z]) {
                stats.funnel_true_in_screen += 1;
            }
            stats.funnel_true_anchor_mates += holds_truth_mates(&anchors);
            // Remembered for the emit site, which is where a half-pair is finally observable but
            // where the candidate list is long gone. A dropped mate either never had an anchor on
            // the true reference (seeding could not offer it) or had one and lost the ranking
            // (selection threw it away) -- and those need opposite fixes.
            truth_anchor_fwd = anchors.iter().any(|AnchorPair(a1, _)| {
                a1.as_ref().is_some_and(|a| a.reference as usize == true_rid)
            });
            truth_anchor_rev = anchors.iter().any(|AnchorPair(_, a2)| {
                a2.as_ref().is_some_and(|a| a.reference as usize == true_rid)
            });
            stats.funnel_true_in_screen_mates += holds_truth_mates(&anchors[0..extend_top_z]);
        }
        let pair_bonus = self.options.args.proper_pair_bonus;
        let (duration, extended_count) = time(|| {
            let extended_count = self.anchor_extender.extend(
                &mut anchors[0..extend_top_z],
                rec_fwd,
                &self.rec_fwd_revc,
                rec_rev,
                &self.rec_rev_revc,
                stats,
            );

            // RESORT before anything downstream selects candidates: extension has just replaced
            // every anchor score with a Hamming-extended one, so the order the extractor produced
            // is stale. The alignment window (top --align-top-y) is taken straight off this, and
            // picking it from a stale order would align the wrong candidates.
            //
            // NB: a reference-id tie-break was tried here and reverted. It made ties deterministic
            // but not more accurate (an arbitrary key cannot favour the truth), and the tuple sort
            // key cost ~60% wall time on the protal marker DB.
            // PROPER-PAIR BONUS. The raw key is score(a1)+score(a2) with a missing side worth 0 --
            // but `core_matches` sums OVERLAPPING seed lengths, so it is not bounded by read length
            // and a seed-rich repeat hit on one mate can outscore a genuine full pair. Measurement:
            // 60% of misclassified pairs and 54,366 half-pairs are a single-sided candidate beating
            // a two-sided one. Two mates agreeing on a locus at the right insert distance is stronger
            // evidence than one mate scoring well, so a pair that HAS both sides is ranked ahead.
            let bonus = pair_bonus;
            glidesort::sort_by_key(&mut anchors[0..extended_count], |AnchorPair(a1, a2)| {
                let s = a1.as_ref().map_or(0, |a| a.score) + a2.as_ref().map_or(0, |a| a.score);
                -(s + if a1.is_some() && a2.is_some() { bonus } else { 0 })
            });

            // assert!(extension_anchors.iter().all(|pair| pair.validate().is_ok()));
            extended_count
        });

        // Assumes sorted anchors !!
        let extension_anchors = &mut anchors[0..extended_count];

        // Mate rescue, LAZY: run only on the candidate(s) that can still be emitted, after extension
        // has re-scored and re-sorted everything.
        //
        // Doing it before extension worked but wasted most of the effort -- 649,934 rescues succeeded
        // to change ~88,000 emitted mates, because a candidate rescued at rank 5 is still rank 5 and
        // never reaches the output. Here the ranking is already final, so a rescue that succeeds is a
        // rescue that counts. The rescued anchor is then extended on its own, so it carries the same
        // extended score and CIGAR as one the index produced and is judged by the same gates.
        if self.options.args.mate_rescue && !extension_anchors.is_empty() {
            let args = &self.options.args;
            let max_insert = 1000u64;
            let max_mm = |len: usize| ((1.0 - args.min_ani) * len as f64) as u32;
            let ext_score = |AnchorPair(a, b): &AnchorPair| {
                a.as_ref().map_or(0, |x| x.score) + b.as_ref().map_or(0, |x| x.score)
            };
            let best_score = ext_score(&extension_anchors[0]);
            // Ties still matter after extension -- selecting rank 1 alone was what made the rank-based
            // trigger insensitive -- but the tail cannot win, so the cap stays small.
            let cutoff = (best_score as f64 * args.mate_rescue_score_frac) as i32;
            let cap = min(args.mate_rescue_top.min(4), extension_anchors.len());
            let full_pair_exists = extension_anchors
                .iter()
                .any(|AnchorPair(a, b)| a.is_some() && b.is_some());
            let mut fired = false;
            let mut rescued_any = false;

            for idx in 0..cap {
                if ext_score(&extension_anchors[idx]) < cutoff {
                    break;
                }
                let AnchorPair(a_fwd, a_rev) = &mut extension_anchors[idx];
                // The partner's orientation follows the ANCHORED mate's strand, not which mate it is.
                // In an FR pair the two mates lie on opposite strands: an anchored mate on the
                // forward strand implies its partner is reverse-complemented, and an anchored mate on
                // the reverse strand implies its partner is forward. Always taking the reverse
                // complement scans the wrong sequence for every read whose anchored mate is reverse
                // -- roughly half of them -- and no offset can then match.
                let (anchored, anchored_len, partner) = match (a_fwd.is_some(), a_rev.is_some()) {
                    (true, false) => {
                        let a = a_fwd.as_ref().unwrap();
                        let p = if a.forward { self.rec_rev_revc.seq() } else { rec_rev.seq() };
                        (a, rec_fwd.seq().len(), p)
                    }
                    (false, true) => {
                        let a = a_rev.as_ref().unwrap();
                        let p = if a.forward { self.rec_fwd_revc.seq() } else { rec_fwd.seq() };
                        (a, rec_rev.seq().len(), p)
                    }
                    _ => continue,
                };
                if anchored.core_matches() < args.mate_rescue_min_core {
                    stats.rescue_skipped_weak += 1;
                    continue;
                }

                stats.rescue_attempted += 1;
                fired = true;
                let mut ambiguous = false;
                let mut indel_suspected = false;
                let rescued = rescue_partner(
                    self.db, anchored, anchored_len, partner, max_insert,
                    max_mm(partner.len()), args.mate_rescue_margin,
                    args.mate_rescue_segments, &mut ambiguous, &mut indel_suspected,
                );
                if ambiguous {
                    stats.rescue_rejected_ambiguous += 1;
                }
                if indel_suspected {
                    stats.rescue_indel_suspected += 1;
                }
                if let Some(a) = rescued {
                    if a_fwd.is_some() { *a_rev = Some(a); } else { *a_fwd = Some(a); }
                    stats.rescue_succeeded += 1;
                    rescued_any = true;
                }
            }
            if fired && full_pair_exists {
                stats.rescue_with_full_pair_available += 1;
            }

            // Extend just what was rescued: a rescued anchor arrives with a position and nothing
            // else, and everything downstream (scores, CIGAR, the emit gates) assumes an extended one.
            if rescued_any {
                let n = min(cap, extension_anchors.len());
                self.anchor_extender.extend(
                    &mut extension_anchors[0..n],
                    rec_fwd,
                    &self.rec_fwd_revc,
                    rec_rev,
                    &self.rec_rev_revc,
                    stats,
                );
            }
        }

        // --debug: after extension + score-sort, is the best anchor the true one? Together with
        // dbg_true_anchor_present this pinpoints the stage: present but not best => extension/scoring
        // ranks a homolog above the truth; not present => the true anchor is never formed (seeding).
        if self.options.args.debug || GOLDSTD_EVAL {
            if true_rid != 0 {
                if let Some(AnchorPair(a1, a2)) = extension_anchors.first() {
                    let best_is_true = a1.as_ref().is_some_and(|a| a.reference as usize == true_rid)
                        || a2.as_ref().is_some_and(|a| a.reference as usize == true_rid);
                    if best_is_true {
                        stats.dbg_true_anchor_best += 1;
                    }
                }
            }
        }

        stats.time_hamming_screen += duration;

        if GOLDSTD_EVAL {
            let mut any_correct = false;
            for AnchorPair(a1, a2) in extension_anchors.iter() {
                let rval = a1.as_ref().map_or_else(|| a2.as_ref().unwrap().reference, |a| a.reference);
                let ref_string = &self.db.get_rname(rval as usize).unwrap();
                let correct = &ref_string.as_bytes()[..min(ref_string.len(), rec_fwd.head().len())]
                    == &rec_fwd.head()[..min(ref_string.len(), rec_fwd.head().len())];
                any_correct |= correct;
            }
            log::trace!("Any correct? {}", any_correct);
        }

        log::trace!(
            "\n##### Sorted and Extended Anchors{}\n\n",
            extension_anchors
                .iter()
                .enumerate()
                .map(|(index, AnchorPair(r1, r2))| {
                    format!(
                        "\n####### Anchor Index: {} \n------Read /1: \n{}\n------Read /2: \n{}",
                        index,
                        if let Some(r1) = r1 {
                            r1.to_string()
                        } else {
                            "No anchor Read /1".to_string()
                        },
                        if let Some(r2) = r2 {
                            r2.to_string()
                        } else {
                            "No anchor Read /2".to_string()
                        },
                    )
                })
                .collect_vec()
                .join("\n")
        );

        //#####################################################################################
        //# At this point enforce valid seeds

        // extension_anchors
        //     .iter_mut()
        //     .enumerate()
        //     .for_each(|(i, (AnchorPair(a1, a2)))| {
        //         match a1 {
        //             Some(a) => {
        //                 if a.seed_status.is_valid()
        //                     && a.seeds.len() > 1
        //                     && a.seeds[0].qbegin() > a.seeds[1].qbegin()
        //                 {
        //                     eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
        //                     panic!("Z 1  {}", a);
        //                 }
        //             }
        //             _ => {}
        //         }

        //         match a2 {
        //             Some(a) => {
        //                 if a.seed_status.is_valid()
        //                     && a.seeds.len() > 1
        //                     && a.seeds[0].qbegin() > a.seeds[1].qbegin()
        //                 {
        //                     eprintln!("{}\n{}\n", rec_fwd.to_string(), rec_rev.to_string());
        //                     panic!("Z 2  {}", a);
        //                 }
        //             }
        //             _ => {}
        //         }
        //     });

        let mut pseudo_mapq = StdPairedAnchorMAPQ::anchor_mapq(extension_anchors);
        let _old_score = &extension_anchors[0].0.as_ref().map_or(0, |a| a.score)
            + &extension_anchors[0].1.as_ref().map_or(0, |a| a.score);


        // Assumes sorted anchors !!
        let anchors_len: usize = extension_anchors.len();
        let mut alignment_anchors =
            &mut extension_anchors[0..min(self.options.args.align_top_y, anchors_len)];

        // Did the true anchor survive the screen's re-ranking into the alignment window? A "no" here
        // with a "yes" at `funnel_in_screen` is the screen mis-ranking the truth -- the failure mode
        // a fixed-diagonal Hamming has on reads with a real indel.
        if GOLDSTD_EVAL && true_rid != 0 {
            if holds_truth(alignment_anchors) {
                stats.funnel_true_in_align_window += 1;
            }
            stats.funnel_true_in_align_window_mates += holds_truth_mates(alignment_anchors);
        }

        // Hoisted out of the closure below, which mutably borrows `self.align`.
        let min_ani = self.options.args.min_ani;

        let (duration, _) = time(|| {
            let mut min_score_1 = None;
            let mut min_score_2 = None;

            // Cost of the best complete PAIR aligned so far, or None until one completes. This is
            // what makes the abort budget dynamic: a candidate that cannot get under it cannot win,
            // so there is no reason to finish aligning it.
            //
            // Bounding each MATE by the pair cost is sound and conservative. Pair score is
            // `s1 + s2` with both <= 0, so beating the best pair needs
            // `s1_new > best_pair - s2_new >= best_pair`, i.e. `cost1_new < best_pair_cost`.
            // A per-mate bound taken from that mate's own best would NOT be sound -- it would
            // discard a weak mate belonging to a winning pair.
            let mut best_pair_cost: Option<i32> = None;
            let mut trace_total = AlignTrace::default();
            let mut stage_hist = [0usize; 6];
            let mut stage_hist_true = [0usize; 6];

            alignment_anchors
                .iter_mut()
                .enumerate()
                .for_each(|(_i, AnchorPair(a1, a2) )| {
                    let reference = match a1 {
                        Some(a) => &self.db.get_reference(a.reference as usize).unwrap(),
                        None => &self
                            .db
                            .get_reference(a2.as_ref().unwrap().reference as usize)
                            .unwrap(),
                    };

                    match a1 {
                        Some(a) => {
                            let query = if a.forward {
                                rec_fwd.seq()
                            } else {
                                self.rec_fwd_revc.seq()
                            };
                            if query.len() == 0 {
                                a.score = 0i32;
                            } else {
                                // Budget the mismatches over the bases that CAN align. Charged over the
                                // whole read, a dovetailing read spends its entire budget on bases outside
                                // the gene and aborts -- measured at 97-98% identity, right marker, a
                                // median 9bp from minibwa's position. Per candidate, not cached: the
                                // alignable window depends on where this anchor sits.
                                let alignable = alignable_len(a, query.len(), reference.len()) as i32;
                                min_score_1 = Some(tighten_budget(
                                    ani_abort_score(min_ani, penalties::MISMATCH, alignable).abs(),
                                    best_pair_cost,
                                    _i,
                                ));
                                self.align.set_max_alignment_score(min_score_1.unwrap());
                                // eprintln!("Align max score: {}", min_score_1.unwrap());

                                if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                                    log::debug!("1  {}", a);
                                }

                                // let status = a.smart_align(&mut self.align, query, reference, 10, min_score_1.unwrap());
                                // let status = a.whole_align(&mut self.align, query, reference, 10, min_score_1.unwrap());

                                let mut trace = AlignTrace::default();
                                let status = a.align(
                                    &mut self.align,
                                    query,
                                    reference,
                                    10,
                                    min_score_1.unwrap(),
                                    screen_hamming(a, query.len()),
                                    &mut trace,
                                );
                                trace_total.merge_from(&trace);
                                if STAGE_TRACE {
                                    stage_hist[trace.stage as usize] += 1;
                                    if GOLDSTD_EVAL && a.reference as usize == true_rid {
                                        stage_hist_true[trace.stage as usize] += 1;
                                    }
                                }

                                // let (qr, rr) = a.whole(query.len(), reference.len());
                                // let (duration, (score, cigar, status)) = time(|| self.align.align(&query[qr], &reference[rr]));

                                match status {
                                    super::common::Status::OK => stats.alignments_successful += 1,
                                    super::common::Status::Dropped => {
                                        stats.alignments_dropped += 1;
                                        // A dropped alignment otherwise keeps its (large, positive)
                                        // pre-alignment extension score, which the score-descending
                                        // re-sort below would then rank ABOVE cleanly-aligned anchors
                                        // (whose scores are <= 0) -- so a divergent homolog that drops
                                        // beats the true locus that aligned. Sink it. (Kept well above
                                        // i32::MIN so the sort's `s1 + s2` cannot overflow.)
                                        a.score = SUNK_SCORE;
                                    }
                                    super::common::Status::Partial => stats.alignments_partial += 1,
                                }

                                let score = a.score;
                                // stats.time_offset += duration;
                                // stats.alignments += 1;
                                // a.score = score / -4;

                                let _ani = 1.0 - a.score as f64 / a.cigar().0.len() as f64 ;
                                // let ani: f64 = (1.0 - a.score as f64/cigar.0.len() as f64);
                                // let ani: f64 = (1.0 - a.score as f64/a.cigar().0.len() as f64);
                                // if score < -50 && score != std::i32::MIN {
                                //     eprintln!("{}/1: {} ANI: {}", i, score, ani);
                                // }

                                // (The tightening that used to live here reassigned min_score_1 from
                                // this candidate's own result, which the next candidate's
                                // unconditional reassignment then discarded -- so it never had any
                                // effect. The cross-candidate bound is now `best_pair_cost`, applied
                                // through `tighten_budget` above.)
                                // eprintln!("{} (asize: {}) Set score {} {} {} {}", i, a.seeds.len(), score, a.score, (1.0 - a.score as f64/cigar.0.len() as f64),  String::from_utf8_lossy(&cigar.0));
                            }
                            // eprintln!("{}", query.len());
                        }
                        None => (),
                    };
                    match a2 {
                        Some(a) => {
                            let query = if a.forward {
                                rec_rev.seq()
                            } else {
                                self.rec_rev_revc.seq()
                            };
                            if query.len() == 0 {
                                a.score = 0i32;
                            } else {
                                // Budget the mismatches over the bases that CAN align. Charged over the
                                // whole read, a dovetailing read spends its entire budget on bases outside
                                // the gene and aborts -- measured at 97-98% identity, right marker, a
                                // median 9bp from minibwa's position. Per candidate, not cached: the
                                // alignable window depends on where this anchor sits.
                                let alignable = alignable_len(a, query.len(), reference.len()) as i32;
                                min_score_2 = Some(tighten_budget(
                                    ani_abort_score(min_ani, penalties::MISMATCH, alignable).abs(),
                                    best_pair_cost,
                                    _i,
                                ));

                                if a.seeds.len() > 1 && a.seeds[0].qbegin() > a.seeds[1].qbegin() {
                                    log::debug!("2  {}", a);
                                }

                                // self.align.set_max_alignment_score(min_score_2.unwrap());
                                // let status = a.smart_align(&mut self.align, query, reference, 10, min_score_2.unwrap());
                                // let status = a.whole_align(&mut self.align, query, reference, 10, min_score_2.unwrap());

                                let mut trace = AlignTrace::default();
                                let status = a.align(
                                    &mut self.align,
                                    query,
                                    reference,
                                    10,
                                    min_score_2.unwrap(),
                                    screen_hamming(a, query.len()),
                                    &mut trace,
                                );
                                trace_total.merge_from(&trace);
                                if STAGE_TRACE {
                                    stage_hist[trace.stage as usize] += 1;
                                    if GOLDSTD_EVAL && a.reference as usize == true_rid {
                                        stage_hist_true[trace.stage as usize] += 1;
                                    }
                                }

                                // let (qr, rr) = a.whole(query.len(), reference.len());
                                // let (duration, (score, cigar, status)) = time(|| self.align.align(&query[qr], &reference[rr]));

                                match status {
                                    super::common::Status::OK => stats.alignments_successful += 1,
                                    super::common::Status::Dropped => {
                                        stats.alignments_dropped += 1;
                                        // A dropped alignment otherwise keeps its (large, positive)
                                        // pre-alignment extension score, which the score-descending
                                        // re-sort below would then rank ABOVE cleanly-aligned anchors
                                        // (whose scores are <= 0) -- so a divergent homolog that drops
                                        // beats the true locus that aligned. Sink it. (Kept well above
                                        // i32::MIN so the sort's `s1 + s2` cannot overflow.)
                                        a.score = SUNK_SCORE;
                                    }
                                    super::common::Status::Partial => stats.alignments_partial += 1,
                                }

                                // match status {
                                //     super::common::Status::OK => {
                                //         if a.reference_cigar_range.len() == 0 {
                                //             eprintln!("Invalid range... {:?}", a.reference_cigar_range);
                                //         }
                                //         if is_alignment_valid(&query, &reference[a.reference_cigar_range.clone()], &a.cigar().0) {
                                //             // print_alignment(&query, &reference[a.reference_cigar_range.clone()], &a.cigar().0);
                                //         } else {
                                //             // eprintln!("------------------------");
                                //             // eprintln!("Valid ? {}", a.validate_seeds(query, reference));
                                //             // eprintln!("{}", a);
                                //             // eprintln!("{}\n{}\n{}",
                                //             //     String::from_utf8_lossy(query),
                                //             //     String::from_utf8_lossy(&reference[a.reference_cigar_range.clone()]),
                                //             //     String::from_utf8_lossy(&a.cigar().0));
                                //             // panic!("Issue {:?}", a.reference_cigar_range.clone());
                                //             eprintln!("Flag issue {:?}", a.reference_cigar_range.clone());
                                //         }
                                //     },
                                //     _ => ()
                                // }

                                let score = a.score;
                                // stats.time_offset += duration;
                                // stats.alignments += 1;
                                // a.score = score / -4;

                                let _ani = 1.0 - a.score as f64 / a.cigar().0.len() as f64 ;
                                // let ani = (1.0 - a.score as f64/cigar.0.len() as f64);
                                // if score < -50 && score != std::i32::MIN {
                                //     eprintln!("{}/2: {} ANI: {}", i, score, ani);
                                // }

                                // (The tightening that used to live here reassigned min_score_2 from
                                // this candidate's own result, which the next candidate's
                                // unconditional reassignment then discarded -- so it never had any
                                // effect. The cross-candidate bound is now `best_pair_cost`, applied
                                // through `tighten_budget` above.)
                                // eprintln!("{} (asize: {}) Set score {} {} {} {}", i, a.seeds.len(), score, a.score, (1.0 - a.score as f64/cigar.0.len() as f64),  String::from_utf8_lossy(&cigar.0));
                            }
                            // eprintln!("{}", query.len());
                        }
                        None => (),
                    };

                    // Both mates of this candidate are done. If the pair completed, its cost becomes
                    // the bound every LATER candidate has to beat (see `best_pair_cost` above).
                    // Sunk anchors (-1_000_000) are excluded: a dropped mate has no real cost.
                    let s1 = a1.as_ref().map_or(0, |a| a.score);
                    let s2 = a2.as_ref().map_or(0, |a| a.score);
                    if s1 > SUNK_SCORE && s2 > SUNK_SCORE {
                        let pair_cost = -(s1 + s2);
                        if pair_cost >= 0 {
                            best_pair_cost =
                                Some(best_pair_cost.map_or(pair_cost, |c| min(c, pair_cost)));
                        }
                    }
                });

            glidesort::sort_by_key(&mut alignment_anchors,|AnchorPair(a1, a2)| {
                let s1 = match a1 {
                    Some(a) => a.score,
                    None => 0,
                };
                let s2 = match a2 {
                    Some(a) => a.score,
                    None => 0,
                };

                - ((s1 + s2) as i64)
            });

            stats.trace.merge_from(&trace_total);
            if STAGE_TRACE {
                for i in 0..stage_hist.len() {
                    stats.align_stage_reached[i] += stage_hist[i];
                    stats.align_stage_reached_true[i] += stage_hist_true[i];
                }
            }
        });
        stats.time_alignment += duration;

        // The anchors have just been re-sorted by ALIGNED score, so best/second-best by alignment
        // are now known -- which is the margin MAPQ should describe. Recompute here rather than
        // keeping the pre-alignment estimate made before any of these candidates were aligned.
        if self.options.args.mapq_from_alignment && !extension_anchors.is_empty() {
            pseudo_mapq = StdPairedAnchorMAPQ::mapq_from(extension_anchors, true);
        }

        if self.options.args.dump_anchors {
            // One write per read: worker threads share stderr, and a single formatted string keeps
            // each read's block intact instead of interleaving line by line.
            let rid = String::from_utf8_lossy(rec_fwd.head()).to_string();
            let rid = rid.split_whitespace().next().unwrap_or("?").trim_start_matches('@')
                .rsplit('/').last().unwrap_or("?").to_string();
            let mut buf = String::new();
            for (rank, AnchorPair(a1, a2)) in extension_anchors.iter().enumerate() {
                for (mate, a) in [(1, a1), (2, a2)] {
                    if let Some(a) = a {
                        let rname = self.db.get_rname(a.reference as usize).unwrap_or("?");
                        buf.push_str(&format!(
                            "ANCHOR\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            rid, mate, rank, rname,
                            a.seeds.first().map_or(0, |sd| sd.rbegin()),
                            a.core_matches() as i32 - a.mismatches as i32,
                            a.score,
                            a.cigar.is_some(),
                        ));
                    }
                }
            }
            if !buf.is_empty() {
                eprint!("{}", buf);
            }
        }

        // eprintln!("SCORE  {} -> {}", old_score, new_score);
        // eprintln!(
        //     "MAPQ   {} -> {}",
        //     pseudo_mapq,
        //     StdPairedAnchorMAPQ::anchor_mapq(extension_anchors)
        // );

        //##########

        let best_extended_anchor_pair = extension_anchors.first().unwrap();

        let reference_id = if best_extended_anchor_pair.0.is_some() {
            &best_extended_anchor_pair.0.as_ref().unwrap().reference
        } else {
            &best_extended_anchor_pair.1.as_ref().unwrap().reference
        };

        let _reference = &self.db.get_reference(*reference_id as usize).unwrap();

        // if pseudo_mapq == 0 && extension_anchors.len() > 1 {
        //     let one = &extension_anchors[0];
        //     let two = &extension_anchors[1];

        //     let rone = if let Some(a) = &one.0 {
        //         a.reference
        //     } else {
        //         one.1.as_ref().unwrap().reference
        //     };
        //     let rtwo = if let Some(a) = &two.0 {
        //         a.reference
        //     } else {
        //         two.1.as_ref().unwrap().reference
        //     };

        //     if rone == rtwo {
        //         eprintln!("Read at fault: {}", String::from_utf8_lossy(rec_fwd.head()));
        //         eprintln!(
        //             "{}\n\nBest Alignment >>\n{}\n{}\n\nSecond Best Alignment>>\n{}{}\n",
        //             pseudo_mapq,
        //             one.0.as_ref().map_or("None".to_string(), |a| a.to_string()),
        //             one.1.as_ref().map_or("None".to_string(), |a| a.to_string()),
        //             two.0.as_ref().map_or("None".to_string(), |a| a.to_string()),
        //             two.1.as_ref().map_or("None".to_string(), |a| a.to_string())
        //         );

        //         if extension_anchors.len() > 2 {
        //             let three = &extension_anchors[2];

        //             let rthree = if let Some(a) = &three.0 {
        //                 a.reference
        //             } else {
        //                 three.1.as_ref().unwrap().reference
        //             };
        //             eprintln!(
        //                 "\nThird Best Alignment >>\n{:?}\n{:?}\n\n",
        //                 three.0, three.1
        //             );
        //         }

        //         let mut input = String::new();
        //         std::io::stdin()
        //             .read_line(&mut input)
        //             .expect("error: unable to read user input");
        //     }

        //     // assert!(rone != rtwo);
        // }

        // ---- Reactivate or put in different piece of code vv

        // let valid_fwd = best_extended_anchor_pair.0.as_ref().map(|a| {
        //     a.validate_seeds(
        //         if a.forward {
        //             rec_fwd.seq()
        //         } else {
        //             self.rec_fwd_revc.seq()
        //         },
        //         reference,
        //     )
        // });

        // log::trace!("MAPQ: {}", pseudo_mapq);

        // let valid_rev = best_extended_anchor_pair.1.as_ref().map(|a| {
        //     a.validate_seeds(
        //         if a.forward {
        //             rec_rev.seq()
        //         } else {
        //             self.rec_rev_revc.seq()
        //         },
        //         reference,
        //     )
        // });
        // let valid = valid_fwd.unwrap_or(true) && valid_rev.unwrap_or(true);

        // ---- Reactivate or put in different piece of code ^^

        // if anchor_pair.0.is_some() {
        //     let a = anchor_pair.0.as_ref().unwrap();
        //     let query = if a.forward { rec_fwd.seq() } else { self.rec_fwd_revc.seq() };
        //     let (qr, rr) = a.whole(query.len(), reference.len());
        //     a.visualize_alignment(query, reference);
        //     eprintln!("Hammingo: {}\n{}\n{}", a.hamming(query, reference), (&query[qr]).ts(),  (&reference[rr]).ts());

        // }
        // if anchor_pair.1.is_some() {
        //     let a = anchor_pair.1.as_ref().unwrap();
        //     let query = if a.forward { rec_rev.seq() } else { self.rec_rev_revc.seq() };
        //     a.visualize_alignment(query, reference);
        //     let (qr, rr) = a.whole(query.len(), reference.len());
        //     eprintln!("Hammingo: {}\n{}\n{}", a.hamming(query, reference), (&query[qr]).ts(),  (&reference[rr]).ts());
        // }

        // if !valid {
        //     eprintln!("Incidence\n{:?} -> {:?}\n{:?} -> {:?}", valid_fwd, anchor_pair.0, valid_rev, anchor_pair.1)
        // }

        // if anchor_pair.0.is_some() && anchor_pair.1.is_some() {
        //     let normal = anchor_pair.0.as_ref().unwrap().forward;
        //     if normal {
        //         eprintln!("{}\n{}", String::from_utf8_lossy(rec_fwd.seq()), String::from_utf8_lossy(self.rec_rev_revc.seq()));
        //     } else {
        //         eprintln!("{}\n{}", String::from_utf8_lossy(self.rec_fwd_revc.seq()), String::from_utf8_lossy(rec_rev.seq()));
        //     }
        // }

        // let before_ref = best_before.reference();
        // let after_ref = best_after.reference();

        // let before_correct = correct(rec_fwd.head(), before_ref, self.db);
        // let after_correct = correct(rec_fwd.head(), after_ref, self.db);

        // if before_correct && !after_correct {
        //     // eprintln!("Anchors... {}", anchors.len());

        //     // eprintln!("Best before: {:?}", best_before);
        //     // eprintln!("Best after:  {:?}", best_after);
        //     anchors.iter().for_each(|AnchorPair(a1, a2)| {
        //         let reference = match a1 {
        //             Some(a) => &self.db.get_reference(a.reference as usize).unwrap(),
        //             None => &self.db.get_reference(a2.as_ref().unwrap().reference as usize).unwrap(),
        //         };

        //         let hamming1 = match a1 {
        //             Some(a) => {
        //                 let query = if a.forward { rec_fwd.seq() } else { self.rec_fwd_revc.seq() };
        //                 if query.len() == 0 { 0 } else {
        //                     query.len() as u64 - a.hamming(query, reference)
        //                 }
        //             },
        //             None => 0,
        //         };
        //         let hamming2 = match a2 {
        //             Some(a) => {
        //                 let query = if a.forward { rec_rev.seq() } else { self.rec_rev_revc.seq() };
        //                 if query.len() == 0 { 0 } else {
        //                     query.len() as u64 - a.hamming(query, reference)
        //                 }
        //             },
        //             None => 0,
        //         };
        //         let score1 = match a1 {
        //             Some(a) => StdAnchorScore::score(a),
        //             None => 0,
        //         };
        //         let score2 = match a2 {
        //             Some(a) => StdAnchorScore::score(a),
        //             None => 0,
        //         };
        //         eprintln!("{}", - ((hamming1 + hamming2) as i64));
        //         eprintln!("anchor {} -- {} {} .....   {} -> {} _____ {} ->{}", AnchorPair(a1.clone(), a2.clone()).reference(),  score2+score1, hamming1+hamming2, score1, hamming1, score2, hamming2);
        //     });

        //     let mut name = String::new();
        //     std::io::stdin().read_line(&mut name).expect("Read line failed.");
        // }

        // ##########################################################################################
        // Output. Both mates are resolved before either is written, because a SAM record carries the
        // other mate's placement, and because --min-ani gates the two mates independently.
        // ##########################################################################################

        let fwd_out = resolve_mate(
            best_extended_anchor_pair.0.as_ref(),
            self.db,
            rec_fwd,
            &self.rec_fwd_revc,
        );
        let rev_out = resolve_mate(
            best_extended_anchor_pair.1.as_ref(),
            self.db,
            rec_rev,
            &self.rec_rev_revc,
        );

        // A mate below the identity OR coverage threshold is treated exactly like a mate with no
        // anchor: nothing is emitted for it, and it does not appear as the other mate's pairing
        // partner either. The coverage floor rejects short local matches to the wrong reference that
        // pass `min_ani` only because the clipped-away remainder of the read costs nothing.
        // End-to-end (clipping only where the read overhangs the reference) replaces the blanket
        // coverage fraction. Setting --min-query-coverage > 0 restores the old gate instead.
        let min_aligned_bases = self.options.args.min_query_coverage;
        // Two independent conditions, both required:
        //   end_to_end     -- WHERE clipping is legitimate: only where the read runs off the
        //                     reference, never in the middle of it. A partial match landing inside
        //                     a gene does not explain the read.
        //   aligned_bases  -- HOW MUCH sequence actually agreed. A read overhanging a marker's end
        //                     by 140 of its 150 bases is legitimately clipped and still carries only
        //                     10 bases of evidence, which is not enough to place it anywhere.
        // Neither alone suffices: a base floor alone accepts mid-reference partial hits, end-to-end
        // alone accepts 10-base edge matches (measured 53.90% precision).
        let require_e2e = !self.options.args.allow_partial;
        let ok = |m: &MateOut| {
            m.ani >= min_ani
                && m.aligned_bases >= min_aligned_bases
                && (!require_e2e || m.end_to_end)
        };
        // PAIR-AWARE base floor. `min_aligned_bases` asks "does this mate carry enough sequence to be
        // placed ON ITS OWN" -- the right question for a single read and the wrong one for a pair. A
        // mate overlapping a marker by 20 bases is unplaceable alone, but if its partner aligns
        // cleanly at the correct insert distance then the pair, jointly, places it: the 20 bases only
        // have to be CONSISTENT with a location the partner already established, not to establish one.
        // Judging the two mates independently discards that evidence and is the same mistake as the
        // pair-level funnel -- treating a pair as two unrelated reads.
        //
        // So a mate that fails only the length test is kept when its partner passed cleanly, down to
        // a lower floor (`--min-query-coverage-mate`). Identity and end-to-end still apply to it in
        // full: this relaxes HOW MUCH sequence is required, never whether the sequence agrees.
        let min_mate_bases = self.options.args.min_query_coverage_mate;
        let ok_backed = |m: &MateOut, partner_ok: bool| {
            partner_ok
                && m.ani >= min_ani
                && m.aligned_bases >= min_mate_bases
                && (!require_e2e || m.end_to_end)
        };
        let fwd_alone = fwd_out.as_ref().is_some_and(&ok);
        let rev_alone = rev_out.as_ref().is_some_and(&ok);
        // Only one side may lean on the other -- otherwise two mates that each fail the full floor
        // could prop each other up on no independent evidence at all.
        let fwd_pass = fwd_alone
            || (min_mate_bases < min_aligned_bases
                && fwd_out.as_ref().is_some_and(|m| ok_backed(m, rev_alone)));
        let rev_pass = rev_alone
            || (min_mate_bases < min_aligned_bases
                && rev_out.as_ref().is_some_and(|m| ok_backed(m, fwd_alone)));
        if fwd_pass && !fwd_alone {
            stats.mate_backed_emitted += 1;
        }
        if rev_pass && !rev_alone {
            stats.mate_backed_emitted += 1;
        }
        if fwd_out.is_some() && !fwd_pass {
            stats.filtered_below_min_ani += 1;
        }
        if rev_out.is_some() && !rev_pass {
            stats.filtered_below_min_ani += 1;
        }
        let fwd_mate = if fwd_pass { fwd_out.as_ref() } else { None };
        let rev_mate = if rev_pass { rev_out.as_ref() } else { None };

        // What SHAPE did this pair emit? Counted always, not only under GOLDSTD_EVAL, because a pair
        // that places one mate and drops the other loses a real alignment while every pair-level
        // counter still calls the pair a success. That asymmetry hid the dominant recall loss:
        // 194k of the 224k mates flexalign missed on markersim had their partner mate emitted.
        match (fwd_mate.is_some(), rev_mate.is_some()) {
            (true, true) => {}
            (false, false) => stats.pairs_none_emitted += 1,
            _ => {
                stats.pairs_half_emitted += 1;
                // Split the half-pair by WHY the dropped mate is missing. Rescue (a scan, not a
                // lookup) can only reach the first case; the second is a ranking failure and needs a
                // different fix entirely, so which one dominates decides where effort goes next.
                if GOLDSTD_EVAL && true_rid != 0 {
                    let dropped_had_anchor = if fwd_mate.is_some() {
                        truth_anchor_rev
                    } else {
                        truth_anchor_fwd
                    };
                    if dropped_had_anchor {
                        stats.halfpair_mate_lost_selection += 1;
                    } else {
                        stats.halfpair_mate_never_anchored += 1;
                    }
                }
            }
        }

        // Bottom of the funnel: the true alignment made it all the way out. Compared against
        // `funnel_in_align_window`, the difference is everything the aligner's own budgets and these
        // two output gates cost -- so a stage that is fine on anchors but lossy on truth shows up as
        // a step down here rather than as a number nobody was counting.
        //
        // Counted BOTH ways. The pair form is satisfied by either mate, so on its own it reports a
        // half-placed pair as fully reported; the mate form (denominator `2 * dbg_checked`) is the
        // one that can see a dropped mate at all.
        if GOLDSTD_EVAL && true_rid != 0 {
            let emitted_truth = |m: &&MateOut| self.db.get_rname(true_rid).is_some_and(|n| n == m.ref_name);
            let fwd_truth = fwd_mate.as_ref().is_some_and(emitted_truth);
            let rev_truth = rev_mate.as_ref().is_some_and(emitted_truth);
            if fwd_truth || rev_truth {
                stats.funnel_true_reported += 1;
            }
            stats.funnel_true_reported_mates += fwd_truth as usize + rev_truth as usize;
        }

        let mut print_reads = false;

        if let Some(m) = fwd_mate {
            if GOLDSTD_EVAL {
                let correct = evaluate::evaluate(
                    stats.gold_std_evaluation.as_mut().unwrap(),
                    m.ref_name,
                    pseudo_mapq as u64,
                    &rec_fwd,
                    self.db,
                    true,
                );
                print_reads |= !correct
            }

            if self.options.args.debug {
                let correct = &m.ref_name.as_bytes()[..min(m.ref_name.len(), rec_fwd.head().len())]
                    == &rec_fwd.head()[..min(m.ref_name.len(), rec_fwd.head().len())];

                if !correct {
                    eprintln!("\n\nIncorrect fwd:");
                    eprintln!("{}", String::from_utf8_lossy(rec_fwd.head()));
                    extension_anchors.print();
                    eprintln!("\nFrom seeds:");
                    eprintln!("\nForward Seeds {}", seeds_fwd.len());
                    for seed in seeds_fwd {
                        eprintln!("\t{}", seed);
                    }
                    eprintln!("\nReverse Seeds {}", seeds_rev.len());
                    for seed in seeds_rev {
                        let seed_ref = self.db.get_rname(seed.rval as usize).unwrap();
                        let seed_correct = &seed_ref.as_bytes()
                            [..min(seed_ref.len(), rec_fwd.head().len())]
                            == &rec_fwd.head()[..min(seed_ref.len(), rec_fwd.head().len())];

                        eprintln!(
                            "\t{} -- {} -- {}",
                            seed,
                            self.db.get_rname(seed.rval as usize).unwrap(),
                            seed_correct
                        );
                    }
                }
            }

            if self.output.has_a() {
                let best = m.anchor;
                self.output.a.as_mut().unwrap().write(
                    &String::from_utf8_lossy(rec_fwd.head()),
                    rec_fwd.seq().len(),
                    best.seeds.first().unwrap().qbegin() as i32,
                    best.seeds.last().unwrap().qend() as i32,
                    best.forward,
                    m.ref_name,
                    m.ref_len,
                    best.seeds.first().unwrap().rbegin() as i32,
                    best.seeds.last().unwrap().rend() as i32,
                    (m.query.len() - m.hamming as usize) as u32,
                    0,
                    pseudo_mapq,
                );
            }

            if self.output.has_b() {
                write_sam_record(
                    self.output.b.as_mut().unwrap(),
                    rec_fwd.head(),
                    m,
                    rev_mate,
                    true,
                    pseudo_mapq,
                    stats,
                );
            }
        } else {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().count_fn();
            }
        }

        if let Some(m) = rev_mate {
            if GOLDSTD_EVAL {
                let correct = evaluate::evaluate(
                    stats.gold_std_evaluation.as_mut().unwrap(),
                    m.ref_name,
                    pseudo_mapq as u64,
                    &rec_rev,
                    self.db,
                    false,
                );
                print_reads |= !correct;
            }

            if self.options.args.debug {
                let correct = &m.ref_name.as_bytes()[..min(m.ref_name.len(), rec_fwd.head().len())]
                    == &rec_fwd.head()[..min(m.ref_name.len(), rec_fwd.head().len())];

                if !correct {
                    eprintln!("\n\nIncorrect Rev:");
                    eprintln!("{}", String::from_utf8_lossy(rec_rev.head()));
                    extension_anchors.print();
                    eprintln!("\nFrom seeds:");
                    eprintln!("\nForward Seeds {}", seeds_fwd.len());
                    for seed in seeds_fwd {
                        let seed_ref = self.db.get_rname(seed.rval as usize).unwrap();
                        let seed_correct = &seed_ref.as_bytes()
                            [..min(seed_ref.len(), rec_rev.head().len())]
                            == &rec_rev.head()[..min(seed_ref.len(), rec_rev.head().len())];

                        eprintln!(
                            "\t{} -- {} -- {}",
                            seed,
                            self.db.get_rname(seed.rval as usize).unwrap(),
                            seed_correct
                        );
                    }
                    eprintln!("\nReverse Seeds {}", seeds_rev.len());
                    for seed in seeds_rev {
                        eprintln!("\t{}", seed);
                    }
                }
            }

            if self.output.has_a() {
                let best = m.anchor;
                self.output.a.as_mut().unwrap().write(
                    &String::from_utf8_lossy(rec_rev.head()),
                    rec_rev.seq().len(),
                    best.seeds.first().unwrap().qbegin() as i32,
                    best.seeds.last().unwrap().qend() as i32,
                    best.forward,
                    m.ref_name,
                    m.ref_len,
                    best.seeds.first().unwrap().rbegin() as i32,
                    best.seeds.last().unwrap().rend() as i32,
                    (m.query.len() - m.hamming as usize) as u32,
                    0,
                    pseudo_mapq,
                );
            }

            if self.output.has_b() {
                write_sam_record(
                    self.output.b.as_mut().unwrap(),
                    rec_rev.head(),
                    m,
                    fwd_mate,
                    false,
                    pseudo_mapq,
                    stats,
                );
            }
        } else {
            if GOLDSTD_EVAL {
                stats.gold_std_evaluation.as_mut().unwrap().count_fn();
            }
        }

        if GOLDSTD_EVAL {
            if print_reads && pseudo_mapq != 0 {
                let eval = stats.gold_std_evaluation.as_mut().unwrap();
                if let Some(or1) = eval.output_read1.as_mut() {
                    assert!(rec_fwd.head().ends_with(b"1"));
                    or1.write(format!("@{}\n", rec_fwd.to_string()));
                }
                if let Some(or2) = eval.output_read2.as_mut() {
                    assert!(rec_rev.head().ends_with(b"2"));
                    or2.write(format!("@{}\n", rec_rev.to_string()));
                }
            }
        }

    }
}


#[cfg(test)]
mod budget_tests {
    use super::*;

    /// The best and the runner-up must never be tightened. MAPQ is the gap between them, so an
    /// early abort on rank 1 truncates the very number that carries the read's ambiguity.
    #[test]
    fn ranks_zero_and_one_keep_the_full_ani_budget() {
        let ani = 300;
        assert_eq!(tighten_budget(ani, Some(10), 0), ani);
        assert_eq!(tighten_budget(ani, Some(10), 1), ani);
    }

    #[test]
    fn later_ranks_are_bounded_by_the_best_pair_cost() {
        let ani = 300;
        assert_eq!(tighten_budget(ani, Some(10), 2), 10, "cannot beat 10, so do not spend 300");
        assert_eq!(tighten_budget(ani, Some(500), 2), ani, "never LOOSEN past --min-ani");
        assert_eq!(tighten_budget(ani, None, 2), ani, "nothing has completed yet");
    }

    /// A perfect pair leaves a zero budget for everyone after it: nothing can beat cost 0, so no
    /// later candidate is worth a single mismatch of work.
    #[test]
    fn a_perfect_pair_closes_the_budget() {
        assert_eq!(tighten_budget(300, Some(0), 2), 0);
    }

    /// `screen_hamming` inverts what the screen stored (`score = query_len - hamming`), and refuses
    /// anything outside [0, query_len] rather than handing the pre-filter a nonsense bound.
    #[test]
    fn screen_hamming_inverts_the_screen_score() {
        let mut a = crate::align::data_structures::anchor::tests::anchor(&[(0, 0, 10)]);
        a.score = 150 - 7;
        assert_eq!(screen_hamming(&a, 150), Some(7));

        a.score = 150;
        assert_eq!(screen_hamming(&a, 150), Some(0), "a perfect screen is zero mismatches");

        // An anchor that never went through the screen carries an unrelated score.
        a.score = -1_000_000;
        assert_eq!(screen_hamming(&a, 150), None, "out of range must not be trusted");
        a.score = 400;
        assert_eq!(screen_hamming(&a, 150), None);
    }
}
