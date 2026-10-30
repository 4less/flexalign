use std::cmp::min;

use bioreader::sequence::fastq_record::{OwnedFastqRecord, RefFastqRecord};

use crate::{
    align::{
        common::{Print},
    },
    database::common::FlexalignDatabase,
    flexalign::stage_time,
    options::Options,
    GOLDSTD_EVAL, STAGE_TRACE,
};

use super::{
    common::{
        penalties, Align, AlignTrace, AnchorExtractor, AnchorPair, AnchorScore, Heuristic,
        KmerExtractor, Or, PAFOutput, PairedAnchorExtender, PairedAnchorExtractor,
        PairedAnchorMAPQ, RangeExtractor, UniqueHit, SAMOutput, SeedExtractor, StdAnchorScore,
        StdPairedAnchorMAPQ,
    },
    data_structures::{anchor::{Anchor, AnchorSeed}, common::Seed},
    process::{
        alignment::{ani_abort_score, approx_ani, hamming_ani},
        evaluate::{self, get_id_from_header},
    },
    sam::{Cigar, Flag},
    stats::Stats,
};
use itertools::Itertools;

use super::process::{alignment::LIBWFA2Alignment, anchor_extender::StdPairedAnchorExtender, anchor_extractor::StdPairedAnchorExtractor};

/// The alignment half of a paired-end pipeline, built one way for every pipeline that aligns: the
/// unsharded run (`process_fastq`) and the sharded join (`shard::disk`, `shard::rejoin_align`).
/// Built separately, the join's had drifted (no anchor gate: other alignments, 4x the reference
/// lookups); built here, an option reaches all of them or none.
pub fn paired_align_half<'a, const K: usize, D: FlexalignDatabase>(
    options: &Options,
    db: &'a D,
) -> (StdPairedAnchorExtractor, StdPairedAnchorExtender<'a, D>, LIBWFA2Alignment) {
    (
        StdPairedAnchorExtractor::new().with_gate(
            (!options.args.no_anchor_gate).then_some(K),
            options.args.anchor_gate_keep,
            options.args.anchor_gate_flank,
        ),
        StdPairedAnchorExtender::new(db),
        LIBWFA2Alignment::default(),
    )
}

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

/// An anchor whose alignment the output can use: present, and not dropped by the aligner.
fn usable(a: &Option<Anchor>) -> Option<&Anchor> {
    a.as_ref().filter(|a| a.score > SUNK_SCORE)
}

/// The best aligned anchor of one mate across the alignment window, for a mate the winning pair
/// has no usable alignment for. Scores are alignment costs (<= 0, higher is better); ties keep
/// the better-ranked pair.
fn best_lone_mate<'p>(window: &'p [AnchorPair], mate: impl Fn(&'p AnchorPair) -> &'p Option<Anchor>) -> Option<&'p Anchor> {
    window
        .iter()
        .filter_map(|p| usable(mate(p)).filter(|a| a.cigar.is_some()))
        .fold(None, |best: Option<&Anchor>, a| match best {
            Some(b) if b.score >= a.score => Some(b),
            _ => Some(a),
        })
}

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

/// Mismatching bytes of `q` against `r` (same length), stopping once they exceed `cap`: eight bytes
/// at a time, the count of non-zero bytes of their XOR. A wrong offset fails within its first word
/// or two, which is what makes trying every start in a fragment's window affordable.
fn mismatches_capped(q: &[u8], r: &[u8], cap: u32) -> u32 {
    const LOW: u64 = 0x0101_0101_0101_0101;
    let mut mm = 0u32;
    let mut qc = q.chunks_exact(8);
    let mut rc = r.chunks_exact(8);
    for (a, b) in (&mut qc).zip(&mut rc) {
        let x = u64::from_ne_bytes(a.try_into().unwrap()) ^ u64::from_ne_bytes(b.try_into().unwrap());
        // A byte is non-zero when any of its bits is: fold them onto its low bit.
        let x = x | (x >> 4);
        let x = x | (x >> 2);
        let x = x | (x >> 1);
        mm += (x & LOW).count_ones();
        if mm > cap {
            return mm;
        }
    }
    mm + qc.remainder().iter().zip(rc.remainder()).filter(|(a, b)| a != b).count() as u32
}

/// The longest fragment a pair is assumed to span, for scans that look for one mate from the other.
const MAX_INSERT: i64 = 1000;

/// OVERHANG RESCUE: finds a mate the index never placed, next to its placed partner, letting it
/// run off the end of the partner's reference.
///
/// On a marker database a fragment often straddles the end of a gene: one mate lies inside it, the
/// other only partly, and the part inside is too short to carry the k-mers seeding needs. The pair
/// geometry still says where it must be -- downstream of a forward partner, upstream of a reverse
/// one, within `MAX_INSERT` -- so every read start in that window is tried, including starts where
/// the read hangs over either end, and each is scored over the bases inside the reference only.
///
/// `mate` is the read as it must map (the reverse complement next to a forward partner). The best
/// start needs at least `min_bases` inside the reference, an identity of `min_ani` there, and a
/// score `margin` above every other start (a mate that fits about equally well at two starts is in
/// a repeat and locates nothing). What it finds is aligned like any candidate, so the identity and
/// base floors downstream judge it as they judge every other mate.
/// The fewest mismatches of `mate` placed gap-free where the pair geometry puts the partner of
/// `placed` (the scan of `rescue_overhang`, without its ambiguity test or the alignment): None when
/// no start is within the identity budget. For comparing pairs, where an ambiguous partner is still
/// a partner.
fn partner_scan_mismatches<D: FlexalignDatabase>(db: &D, placed: &Anchor, mate: &[u8], min_ani: f64) -> Option<u32> {
    let reference = db.get_reference(placed.reference as usize)?;
    let (n, len) = (reference.len() as i64, mate.len() as i64);
    if len == 0 || placed.reference_cigar_range.is_empty() {
        return None;
    }
    let (ps, pe) = (placed.reference_cigar_range.start as i64, placed.reference_cigar_range.end as i64);
    let (lo, hi) = if placed.forward { (ps, ps + MAX_INSERT - len) } else { (pe - MAX_INSERT, pe - len) };
    let (lo, hi) = (lo.max(0), hi.min(n - len));
    let budget = ((1.0 - min_ani) * len as f64) as u32;
    let mut best: Option<u32> = None;
    for st in lo..=hi {
        let r = &reference[st as usize..(st + len) as usize];
        let mm = mismatches_capped(mate, r, best.map_or(budget, |b| b.min(budget)));
        if mm <= budget && best.is_none_or(|b| mm < b) {
            best = Some(mm);
            if mm == 0 {
                break;
            }
        }
    }
    best
}

fn rescue_overhang<D: FlexalignDatabase>(
    db: &D,
    aligner: &mut (impl Align + Heuristic),
    placed: &Anchor,
    mate: &[u8],
    min_ani: f64,
    min_bases: usize,
    margin: i32,
) -> Option<Anchor> {
    let reference = db.get_reference(placed.reference as usize)?;
    let (n, len) = (reference.len() as i64, mate.len() as i64);
    if len == 0 || placed.reference_cigar_range.is_empty() || min_bases == 0 {
        return None;
    }
    let (ps, pe) = (placed.reference_cigar_range.start as i64, placed.reference_cigar_range.end as i64);
    let min_bases = min_bases as i64;
    // Read starts the geometry allows, clipped to those leaving `min_bases` inside the reference.
    let (lo, hi) = if placed.forward { (ps, ps + MAX_INSERT - len) } else { (pe - MAX_INSERT, pe - len) };
    let (lo, hi) = (lo.max(min_bases - len), hi.min(n - min_bases));
    if lo > hi {
        return None;
    }
    // Score of a start: matches - 4 x mismatches over the overlap, bwa's match/mismatch ratio.
    let mut best: Option<(i64, i64, u32, i64)> = None; // (score, start, mismatches, overlap)
    let mut second = i64::MIN;
    for st in lo..=hi {
        let (rs, re) = (st.max(0), (st + len).min(n));
        let ov = re - rs;
        let budget = ((1.0 - min_ani) * ov as f64) as u32;
        let q = &mate[(rs - st) as usize..(re - st) as usize];
        let r = &reference[rs as usize..re as usize];
        let mm = mismatches_capped(q, r, budget);
        if mm > budget {
            continue;
        }
        let score = ov - 5 * mm as i64;
        match best {
            Some((b, ..)) if score <= b => second = second.max(score),
            _ => {
                if let Some((b, ..)) = best {
                    second = second.max(b);
                }
                best = Some((score, st, mm, ov));
            }
        }
    }
    let (score, st, mm, ov) = best?;
    if second != i64::MIN && score - second < margin as i64 {
        return None;
    }
    // Seeds are exact matches -- the aligner takes what lies between two seeds as the part to align
    // -- so the overlap goes in as its runs of agreement (those of 8+ bases; the longest run when
    // none is), not as one seed with mismatches inside it. A reverse anchor takes its seeds in
    // FORWARD-read coordinates and flips them itself (`set_forward`), so for one the offsets are
    // handed over mirrored.
    let (rs, qoff) = (st.max(0), st.max(0) - st);
    let forward = !placed.forward;
    let q = &mate[qoff as usize..(qoff + ov) as usize];
    let r = &reference[rs as usize..(rs + ov) as usize];
    let mut runs: Vec<(i64, i64)> = Vec::new(); // (offset in the overlap, length)
    let mut i = 0i64;
    while i < ov {
        if q[i as usize] != r[i as usize] {
            i += 1;
            continue;
        }
        let start = i;
        while i < ov && q[i as usize] == r[i as usize] {
            i += 1;
        }
        runs.push((start, i - start));
    }
    let longest = runs.iter().map(|r| r.1).max()?;
    let seeds: Vec<(u64, i64, i64)> = runs
        .iter()
        .filter(|&&(_, l)| l >= 8.min(longest))
        .map(|&(o, l)| {
            let qpos = qoff + o;
            ((rs + o) as u64, if forward { qpos } else { len - qpos - l }, l)
        })
        .collect();
    let (r0, q0, l0) = seeds[0];
    let mut a = Anchor::from_seed(&Seed {
        rpos: r0,
        rval: placed.reference,
        qpos: q0 as u32,
        mismatch: 0,
        length: l0.min(u8::MAX as i64) as u8,
        flag: 0,
    });
    for &(rp, qp, l) in &seeds[1..] {
        a.seeds.push(AnchorSeed { qpos: qp as u32, rpos: rp, length: l as u32 });
        a.seed_count += 1;
    }
    a.mismatches = mm;
    // Not evidence from the index.
    a.counter1 = a.seed_count as u16;
    a.set_forward(forward, mate.len());
    // What the extender establishes for an anchor from the index (`anchor_extender`): seeds in
    // order on both sequences, which the seed-based alignment requires. Runs of one diagonal are
    // ordered by construction.
    debug_assert!(a.seeds.windows(2).all(|w| w[0].qend() <= w[1].qbegin() && w[0].rend() <= w[1].rbegin()));
    a.seed_status.set_valid();
    let alignable = alignable_len(&a, mate.len(), reference.len()) as i32;
    let budget = ani_abort_score(min_ani, penalties::MISMATCH, alignable).abs();
    match a.align(aligner, mate, reference, 10, budget, None, &mut AlignTrace::default()) {
        super::common::Status::Dropped => None,
        _ => Some(a),
    }
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

/// Buffers of the fast-hit path, reused across pairs.
#[derive(Clone, Default)]
pub struct FastHitScratch {
    hits_fwd: Vec<UniqueHit>,
    hits_rev: Vec<UniqueHit>,
    cigar_fwd: Vec<u8>,
    cigar_rev: Vec<u8>,
    votes: Vec<(u64, i64, bool, u32)>,
    /// Open-addressing index into `votes` (index + 1; 0 is empty), reused across mates.
    vote_slots: Vec<u32>,
    cands_fwd: Vec<FastCand>,
    cands_rev: Vec<FastCand>,
    /// (key, index) table of the anchor-pair sorts, reused across reads.
    sort_scratch: Vec<(i64, u32)>,
}

/// One candidate placement of a mate on the fast-hit path: the read lies gapless on
/// `[start, end)` of `rval`, `mm` mismatches, `votes` k-mers on that diagonal, and `score` is
/// `K * votes - mm`, the anchor score (`core_matches - mismatches`) the full path ranks by.
#[derive(Clone, Copy)]
struct FastCand {
    rval: u64,
    start: usize,
    end: usize,
    forward: bool,
    mm: usize,
    score: i64,
}

/// Why a mate fell through to the full path.
enum FastMiss {
    FewVotes,
    OffReference,
    Mismatches,
}

/// Lists a mate's candidate placements from its few-position k-mers, or says why there are none.
///
/// `K` is the k-mer length (the extractor's `C` parameter in `ModularPE`). Each hit votes for a
/// (reference, diagonal) under both strand hypotheses -- read forward (`rpos - qpos`) and read
/// reverse-complement (`rpos - (len - K - qpos)`); the index stores the canonical k-mer, so a hit
/// says nothing about strand by itself. The true strand's votes pile onto one diagonal while the
/// other strand's scatter. Every locus with at least `min_hits` votes and at least half the top
/// count is a candidate: the oriented read is compared base by base against its window (which
/// rejects a wrong strand call outright); candidates inside the reference and within
/// `max_mismatch` are returned, best score first. The pairing then picks between them, as the
/// full path pairs anchors, so two copies of a locus are told apart by mismatches AND the mate.
#[allow(clippy::too_many_arguments)]
fn fast_hit_mate<const K: usize, D: FlexalignDatabase>(
    db: &D,
    hits: &[UniqueHit],
    rec: &RefFastqRecord,
    revc: &mut OwnedFastqRecord,
    votes: &mut Vec<(u64, i64, bool, u32)>,
    slots: &mut Vec<u32>,
    cands: &mut Vec<FastCand>,
    min_hits: usize,
    max_mismatch: usize,
) -> Result<(), FastMiss> {
    cands.clear();
    let len = rec.seq().len();
    if hits.len() < min_hits || len < K {
        return Err(FastMiss::FewVotes);
    }
    votes.clear();
    // Votes in first-seen order (the candidates' tie order), counted through a hash table of at
    // most half load: two votes a hit, four slots a hit.
    let size = (4 * hits.len()).next_power_of_two();
    if slots.len() < size {
        slots.resize(size, 0);
    }
    let slots = &mut slots[..size];
    slots.fill(0);
    let mask = size - 1;
    for h in hits {
        let fwd = (h.rval, h.rpos as i64 - h.qpos as i64, true);
        let rev = (h.rval, h.rpos as i64 - (len - K - h.qpos as usize) as i64, false);
        for key in [fwd, rev] {
            let hash = (key.0 ^ (key.1 as u64).rotate_left(24) ^ (key.2 as u64) << 63).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            let mut s = (hash >> 32) as usize & mask;
            loop {
                match slots[s] {
                    0 => {
                        slots[s] = votes.len() as u32 + 1;
                        votes.push((key.0, key.1, key.2, 1));
                        break;
                    }
                    v => {
                        let e = &mut votes[v as usize - 1];
                        if (e.0, e.1, e.2) == key {
                            e.3 += 1;
                            break;
                        }
                    }
                }
                s = (s + 1) & mask;
            }
        }
    }
    let top = votes.iter().map(|v| v.3).max().unwrap_or(0);
    if (top as usize) < min_hits {
        return Err(FastMiss::FewVotes);
    }
    // Touch every candidate's window first so the reference misses overlap, then evaluate.
    for &(rval, diag, _, n) in votes.iter() {
        if (n as usize) < min_hits || 2 * n < top || diag < 0 {
            continue;
        }
        if let Some(reference) = db.get_reference(rval as usize) {
            let start = diag as usize;
            if start + len <= reference.len() {
                let w = &reference[start..start + len];
                let mut off = 0;
                while off < w.len() {
                    prefetch_line(&w[off..]);
                    off += 64;
                }
            }
        }
    }
    let mut revc_ready = false;
    let mut off_reference = 0usize;
    let mut over_cap = 0usize;
    for &(rval, diag, forward, n) in votes.iter() {
        if (n as usize) < min_hits || 2 * n < top {
            continue;
        }
        let reference = match db.get_reference(rval as usize) {
            Some(r) => r,
            None => continue,
        };
        if diag < 0 || diag as usize + len > reference.len() {
            off_reference += 1;
            continue;
        }
        let start = diag as usize;
        let window = &reference[start..start + len];
        let query: &[u8] = if forward {
            rec.seq()
        } else {
            if !revc_ready {
                rec.reverse_complement(revc);
                revc_ready = true;
            }
            revc.seq()
        };
        let mm = triple_accel::hamming(query, window) as usize;
        if mm > max_mismatch {
            over_cap += 1;
            continue;
        }
        cands.push(FastCand { rval, start, end: start + len, forward, mm, score: K as i64 * n as i64 - mm as i64 });
    }
    if cands.is_empty() {
        return Err(if over_cap > 0 { FastMiss::Mismatches } else if off_reference > 0 { FastMiss::OffReference } else { FastMiss::FewVotes });
    }
    cands.sort_by_key(|c| (c.mm, -c.score));
    Ok(())
}

/// Whether two placements make a pair: one reference, opposite strands, within `max_insert`.
fn fast_pairable(a: &FastCand, b: &FastCand, max_insert: usize) -> bool {
    a.rval == b.rval && a.forward != b.forward && a.end.max(b.end) - a.start.min(b.start) <= max_insert
}

/// The fast-hit path for one pair (`--fast-hit`): both mates' candidates from `fast_hit_mate`,
/// then the pairable combination with the fewest mismatches (most votes on a tie), written as a
/// proper pair. MAPQ is the pair's score margin over the runner-up -- the best other pairable
/// combination, or the best lone candidate on either mate -- scaled as `StdPairedAnchorMAPQ`
/// scales the anchor margin (60 for no runner-up). True when the pair was written; false leaves
/// it to the full path.
#[allow(clippy::too_many_arguments)]
fn fast_hit_pair<const K: usize, const F: usize, RE: RangeExtractor<K, F>, SO: SAMOutput, D: FlexalignDatabase>(
    args: &crate::options::Args,
    db: &D,
    re_fwd: &mut RE,
    re_rev: &mut RE,
    kmers_fwd: &[(usize, kmerrs::consecutive::kmer::Kmer<K>)],
    kmers_rev: &[(usize, kmerrs::consecutive::kmer::Kmer<K>)],
    rec_fwd: &RefFastqRecord,
    rec_rev: &RefFastqRecord,
    revc_fwd: &mut OwnedFastqRecord,
    revc_rev: &mut OwnedFastqRecord,
    scratch: &mut FastHitScratch,
    out: &mut SO,
    stats: &mut Stats,
    write: bool,
) -> bool {
    const MAX_INSERT: usize = 1000;
    let (min_hits, max_mm, max_pos) = (args.fast_hit_min_hits, args.fast_hit_max_mismatch, args.fast_hit_max_positions);
    let stride = args.fast_hit_stride.max(1);
    let mut miss = |m: FastMiss, stats: &mut Stats| {
        match m {
            FastMiss::FewVotes => stats.fast_hit_few_votes += 1,
            FastMiss::OffReference => stats.fast_hit_off_reference += 1,
            FastMiss::Mismatches => stats.fast_hit_mismatches += 1,
        }
        false
    };
    // Two rounds: every `stride`-th k-mer first, which places most pairs; the full set when that
    // leaves a mate short of votes (any other miss is final: more k-mers change neither the
    // window nor its mismatches). The full round's ranges are kept for the fallback's generate().
    let mut listed = false;
    for round in 0..2 {
        let full = round == 1 || stride == 1;
        let ok = if full {
            re_fwd.unique_hits(kmers_fwd, max_pos, &mut scratch.hits_fwd)
        } else {
            re_fwd.unique_hits_strided(kmers_fwd, stride, max_pos, &mut scratch.hits_fwd)
        };
        if !ok {
            return false;
        }
        match fast_hit_mate::<K, D>(db, &scratch.hits_fwd, rec_fwd, revc_fwd, &mut scratch.votes, &mut scratch.vote_slots, &mut scratch.cands_fwd, min_hits, max_mm) {
            Ok(()) => {}
            Err(FastMiss::FewVotes) if !full => continue,
            Err(e) => { if !full { re_fwd.unique_hits(kmers_fwd, max_pos, &mut scratch.hits_fwd); } return miss(e, stats); }
        }
        let ok = if full {
            re_rev.unique_hits(kmers_rev, max_pos, &mut scratch.hits_rev)
        } else {
            re_rev.unique_hits_strided(kmers_rev, stride, max_pos, &mut scratch.hits_rev)
        };
        if !ok {
            return false;
        }
        match fast_hit_mate::<K, D>(db, &scratch.hits_rev, rec_rev, revc_rev, &mut scratch.votes, &mut scratch.vote_slots, &mut scratch.cands_rev, min_hits, max_mm) {
            Ok(()) => {}
            Err(FastMiss::FewVotes) if !full => continue,
            Err(e) => { if !full { re_fwd.unique_hits(kmers_fwd, max_pos, &mut scratch.hits_fwd); re_rev.unique_hits(kmers_rev, max_pos, &mut scratch.hits_rev); } return miss(e, stats); }
        }
        listed = true;
        break;
    }
    if !listed {
        return miss(FastMiss::FewVotes, stats);
    }
    // The pairable combination with the fewest mismatches; most votes (score) on a tie.
    let (c1s, c2s) = (&scratch.cands_fwd, &scratch.cands_rev);
    let mut best: Option<(usize, usize, usize, i64)> = None; // (i, j, mm, score)
    let mut runner_up: i64 = 0;
    for (i, a) in c1s.iter().enumerate() {
        for (j, b) in c2s.iter().enumerate() {
            if !fast_pairable(a, b, MAX_INSERT) {
                continue;
            }
            let (mm, score) = (a.mm + b.mm, a.score + b.score);
            match best {
                Some((_, _, bmm, bscore)) if (mm, -score) >= (bmm, -bscore) => runner_up = runner_up.max(score),
                Some((_, _, _, bscore)) => { runner_up = runner_up.max(bscore); best = Some((i, j, mm, score)); }
                None => best = Some((i, j, mm, score)),
            }
        }
    }
    let (i, j, _, score) = match best {
        Some(b) => b,
        None => { stats.fast_hit_unpaired += 1; return false; }
    };
    let (m1, m2) = (c1s[i], c2s[j]);
    if !write {
        return true; // the probe only asks whether the pair would have been placed
    }
    // A lone candidate elsewhere competes too, as a single-sided anchor pair would.
    for (k, c) in c1s.iter().enumerate() { if k != i { runner_up = runner_up.max(c.score); } }
    for (k, c) in c2s.iter().enumerate() { if k != j { runner_up = runner_up.max(c.score); } }
    let ref_name = match db.get_rname(m1.rval as usize) {
        Some(n) => n,
        None => return false,
    };
    // As StdPairedAnchorMAPQ::mapq_from scores anchor pairs: margin of the pair's score over the
    // runner-up's, scaled to 60.
    let s1 = score.max(1) as f64;
    let s2 = runner_up.max(0) as f64;
    let mut mapq = (60.0 * (1.0 - s2 / s1)).round().clamp(0.0, 60.0) as u8;
    // With --mapq-per-mismatch: the margin in mismatches over the runner-up -- another pairable
    // combination, or a lone candidate charged the unpaired penalty -- as the full path takes it
    // from the alignments. Votes tie whenever a core has its two positions on two copies; the
    // mismatches at the two copies are what tell them apart.
    let per_mm = crate::align::common::MAPQ_PER_MISMATCH.load(std::sync::atomic::Ordering::Relaxed) as i64;
    if per_mm > 0 {
        let mmc = penalties::MISMATCH as i64;
        let unpaired = crate::align::common::UNPAIRED_PENALTY.load(std::sync::atomic::Ordering::Relaxed) as i64;
        let best_cost = (m1.mm + m2.mm) as i64 * mmc;
        let mut alt = i64::MAX;
        for (a_i, a) in c1s.iter().enumerate() {
            for (b_i, b) in c2s.iter().enumerate() {
                if (a_i, b_i) != (i, j) && fast_pairable(a, b, MAX_INSERT) {
                    alt = alt.min((a.mm + b.mm) as i64 * mmc);
                }
            }
        }
        for (k, c) in c1s.iter().enumerate() {
            if k != i { alt = alt.min(c.mm as i64 * mmc + m2.mm as i64 * mmc + unpaired); }
        }
        for (k, c) in c2s.iter().enumerate() {
            if k != j { alt = alt.min(c.mm as i64 * mmc + m1.mm as i64 * mmc + unpaired); }
        }
        if alt != i64::MAX {
            mapq = (per_mm * ((alt - best_cost).max(0) / mmc)).min(60) as u8;
        }
    }
    if dump_this(rec_fwd.head()) {
        // --dump-candidates: an F line per fast-placed pair. alt_pair: the fewest mismatches of
        // another pairable combination; alt_lone: of a lone candidate plus the placed partner's.
        let mut alt_pair = usize::MAX;
        for (a_i, a) in c1s.iter().enumerate() {
            for (b_i, b) in c2s.iter().enumerate() {
                if (a_i, b_i) != (i, j) && fast_pairable(a, b, MAX_INSERT) {
                    alt_pair = alt_pair.min(a.mm + b.mm);
                }
            }
        }
        let mut alt_lone = usize::MAX;
        for (k, c) in c1s.iter().enumerate() { if k != i { alt_lone = alt_lone.min(c.mm + m2.mm); } }
        for (k, c) in c2s.iter().enumerate() { if k != j { alt_lone = alt_lone.min(c.mm + m1.mm); } }
        let opt = |x: usize| if x == usize::MAX { "-".to_string() } else { x.to_string() };
        let name = String::from_utf8_lossy(rec_fwd.head()).split_whitespace().next().unwrap_or("").trim_start_matches('@').to_string();
        let name = name.strip_suffix("/1").unwrap_or(&name).to_string();
        let line = format!("F\t{name}\t{mapq}\t{ref_name}:{}\t{ref_name}:{}\t{}\t{}\t{}\t{}\t{score}\t{runner_up}\t{}\t{}\n",
            m1.start, m2.start, m1.mm, m2.mm, c1s.len(), c2s.len(), opt(alt_pair), opt(alt_lone));
        if let Some((w, _)) = DUMP.get() {
            use std::io::Write as _;
            if let Some(w) = w.lock().unwrap().as_mut() {
                let _ = w.write_all(line.as_bytes());
            }
        }
    }
    let leftmost = m1.start.min(m2.start);
    let rightmost = m1.end.max(m2.end);
    let span = (rightmost - leftmost) as i64;
    let revc_fwd: &OwnedFastqRecord = revc_fwd;
    let revc_rev: &OwnedFastqRecord = revc_rev;
    let mut cigar_fwd = Cigar(std::mem::take(&mut scratch.cigar_fwd));
    let mut cigar_rev = Cigar(std::mem::take(&mut scratch.cigar_rev));
    for (first, m, mate, cigar, rec, revc) in [
        (true, &m1, &m2, &mut cigar_fwd, rec_fwd, revc_fwd),
        (false, &m2, &m1, &mut cigar_rev, rec_rev, revc_rev),
    ] {
        let (seq, qual) = if m.forward { (rec.seq(), rec.qual()) } else { (revc.seq(), revc.qual()) };
        let reference = match db.get_reference(m.rval as usize) { Some(r) => r, None => return false };
        let window = &reference[m.start..m.end];
        cigar.0.clear();
        cigar.0.extend(seq.iter().zip(window).map(|(q, r)| if q == r { b'M' } else { b'X' }));
        let mut flag = Flag::new();
        flag.paired_end(true);
        flag.proper_pair(true);
        if first { flag.first_in_pair(true); } else { flag.last_in_pair(true); }
        if !m.forward { flag.reverse(true); }
        if !mate.forward { flag.mate_reverse(true); }
        // The leftmost mate reports the span as positive, the other as negative.
        let tlen = if m.start <= mate.start { span } else { -span };
        out.write(sam_qname(rec.head()), flag.value(), ref_name, m.start, mapq, cigar, Some(ref_name), mate.start, tlen, seq, qual);
    }
    scratch.cigar_fwd = cigar_fwd.0;
    scratch.cigar_rev = cigar_rev.0;
    true
}

/// Writes one SAM record for `m`, using `mate` for the pairing columns.
/// `--discordant-mapq-cap`, set once at start-up (255: off).
pub static DISCORDANT_MAPQ_CAP: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(255);
/// `--discordant-only-mapped`: the cap only for a mate whose partner is written elsewhere.
pub static DISCORDANT_ONLY_MAPPED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn write_sam_record<SO: SAMOutput>(
    out: &mut SO,
    head: &[u8],
    m: &MateOut,
    mate: Option<&MateOut>,
    first_in_pair: bool,
    mapping_quality: u8,
    lone_mate_penalty: u8,
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

    // `--lone-mate-mapq-penalty`: a mate without its partner on the same reference. It stays at 1 or
    // above when it had MAPQ, so it still ranks over the calls that were ambiguous to begin with.
    let supported = mate.is_some_and(|mm| !mm.anchor.reference_cigar_range.is_empty() && mm.anchor.reference == anchor.reference);
    let mapping_quality = if supported || mapping_quality == 0 { mapping_quality } else { mapping_quality.saturating_sub(lone_mate_penalty).max(1) };

    // SAM reserves MAPQ 255 for "unavailable". The pipeline's pseudo-mapq is a raw score gap, not a
    // phred-scaled probability, and can land on 255 -- clamp so it is never mistaken for "missing".
    let mapping_quality = mapping_quality.min(254);
    // `--discordant-mapq-cap`: a mate whose partner is unmapped or not within MAX_INSERT of it. On a
    // genome these are mostly wrong at any MAPQ (human: 7% of them right at MAPQ >= 1).
    let together = mate.is_some_and(|mm| {
        !mm.anchor.reference_cigar_range.is_empty()
            && mm.anchor.reference == anchor.reference
            && mm.anchor.reference_cigar_range.start.abs_diff(anchor.reference_cigar_range.start) < MAX_INSERT as usize
    });
    let cap = DISCORDANT_MAPQ_CAP.load(std::sync::atomic::Ordering::Relaxed);
    // A partner not written at all may lie off the reference's end (a marker gene's): not evidence
    // against this mate with --discordant-only-mapped.
    let partner_elsewhere = mate.is_some_and(|mm| !mm.anchor.reference_cigar_range.is_empty());
    let exempt = DISCORDANT_ONLY_MAPPED.load(std::sync::atomic::Ordering::Relaxed) && !partner_elsewhere;
    let mapping_quality = if together || exempt { mapping_quality } else { mapping_quality.min(cap) };

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
        // A read (pair) starts: no reference slice of the previous one is used again.
        crate::database::refseq2::reset_scratch();
        stats.reads_processed += 1;

        let (duration, kmers) = stage_time(|| self.kmer_extractor.generate(rec, stats));
        stats.time_get_kmers += duration;

        let (duration, ranges) = stage_time(|| self.range_extractor.generate(kmers, rec, stats));
        stats.time_get_ranges += duration;

        let (duration, seeds) = stage_time(|| self.seed_extractor.generate(ranges, stats));
        stats.time_range_header += duration;
        stats.seeds += seeds.len();

        let (duration, anchors) = stage_time(|| {
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

        let (_duration, _) = stage_time(|| {
            anchors.sort_unstable_by_key(|a| {
                -((a.core_matches() - a.mismatches as usize - a.indels()) as i64)
            });
        });

        let (duration, _) = stage_time(|| {
            rec.reverse_complement(&mut self.rec_rev);
        });
        stats.time_reverse_complement += duration;
        stats.time_anchor_sorting += duration;

        let best = anchors.first().unwrap();
        let ref_string = &self.db.get_rname(best.reference as usize).unwrap();
        let reference = self.db.get_reference(best.reference as usize).unwrap();

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
    /// Scratch of the fast-hit path (`--fast-hit`), reused across pairs.
    pub fast: FastHitScratch,
}

/// How many pairs the probe tries the fast-hit path on before judging it, and the success rate
/// (in fifths) it must reach to stay on. On a marker-gene reference only a few % of pairs have
/// any home, so nearly every attempt falls through and the attempt itself is pure overhead
/// (~5% of a run); on a whole-genome reference three quarters succeed. The probe runs once, on
/// the first pairs of the input, before the workers start -- a per-worker probe depended on
/// which batches each thread happened to get and made the output differ between runs.
pub const FAST_HIT_PROBE: usize = 2048;
pub const FAST_HIT_MIN_FIFTHS: usize = 1;

/// The probe's verdict for this run; `--no-fast-hit` is checked separately.
pub static FAST_HIT_ENABLED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(true);

pub fn set_fast_hit_enabled(on: bool) {
    FAST_HIT_ENABLED.store(on, std::sync::atomic::Ordering::Relaxed);
}

pub fn fast_hit_enabled() -> bool {
    FAST_HIT_ENABLED.load(std::sync::atomic::Ordering::Relaxed)
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
    /// The probe: would the fast-hit path place this pair? Writes nothing.
    pub fn probe_fast_hit(&mut self, rec_fwd: &RefFastqRecord, rec_rev: &RefFastqRecord, stats: &mut Stats) -> bool {
        // A read (pair) starts: no reference slice of the previous one is used again.
        crate::database::refseq2::reset_scratch();
        if !self.output.has_b() {
            return false;
        }
        let kmers_fwd = self.kmer_extractor_fwd.generate(rec_fwd, stats);
        let kmers_rev = self.kmer_extractor_rev.generate(rec_rev, stats);
        fast_hit_pair::<C, F, RE, SO, D>(
            &self.options.args,
            self.db,
            &mut self.range_extractor_fwd,
            &mut self.range_extractor_rev,
            kmers_fwd,
            kmers_rev,
            rec_fwd,
            rec_rev,
            &mut self.rec_fwd_revc,
            &mut self.rec_rev_revc,
            &mut self.fast,
            self.output.b.as_mut().unwrap(),
            stats,
            false,
        )
    }

    pub fn run(
        &mut self,
        rec_fwd: &RefFastqRecord,
        rec_rev: &RefFastqRecord,
        stats: &mut Stats,
    ) -> () {
        if STAGE_TRACE {
            let started = std::time::Instant::now();
            self.run_pair(rec_fwd, rec_rev, stats);
            stats.trace_counts.time_pair += started.elapsed();
        } else {
            self.run_pair(rec_fwd, rec_rev, stats);
        }
    }

    fn run_pair(
        &mut self,
        rec_fwd: &RefFastqRecord,
        rec_rev: &RefFastqRecord,
        stats: &mut Stats,
    ) -> () {
        // A read (pair) starts: no reference slice of the previous one is used again.
        crate::database::refseq2::reset_scratch();
        stats.reads_processed += 2;

        log::trace!("\n####################################################################\n#### New pair \n####################################################################");

        if GOLDSTD_EVAL {
            let header_str = String::from_utf8_lossy(rec_fwd.head());
            let id = get_id_from_header(&header_str, self.db);
            log::trace!("True ID: {}\n", id);
        }

        log::trace!("{}\n{}", rec_fwd.to_string(), rec_rev.to_string());

        // Extract minimizer. Each mate's key lines are fetched as soon as its k-mers are known, so
        // the second mate's extraction hides the first's misses: one read's keys are too few for
        // the lookup's own fetch-ahead to cover them.
        let (duration, kmers_fwd) = stage_time(|| self.kmer_extractor_fwd.generate(rec_fwd, stats));
        let n_kmers_fwd = kmers_fwd.len();
        stats.time_get_kmers += duration;
        self.range_extractor_fwd.prefetch(kmers_fwd);
        let (duration, kmers_rev) = stage_time(|| self.kmer_extractor_rev.generate(rec_rev, stats));
        let n_kmers_rev = kmers_rev.len();
        stats.time_get_kmers += duration;
        self.range_extractor_rev.prefetch(kmers_rev);

        // The fast-hit path: a pair pinned by its unique k-mers is written here and never seeded.
        if self.options.args.fast_hit_on() && fast_hit_enabled() && self.output.has_b() {
            let fast_started = STAGE_TRACE.then(std::time::Instant::now);
            let done = fast_hit_pair::<C, F, RE, SO, D>(
                &self.options.args,
                self.db,
                &mut self.range_extractor_fwd,
                &mut self.range_extractor_rev,
                kmers_fwd,
                kmers_rev,
                rec_fwd,
                rec_rev,
                &mut self.rec_fwd_revc,
                &mut self.rec_rev_revc,
                &mut self.fast,
                self.output.b.as_mut().unwrap(),
                stats,
                true,
            );
            if let Some(t) = fast_started {
                stats.trace_counts.fast_hit_tried += 1;
                stats.trace_counts.time_fast_hit += t.elapsed();
            }
            if done {
                stats.fast_hit_pairs += 1;
                return;
            }
        }

        // Get ranges from minimizers
        let (duration, ranges_fwd) =
            stage_time(|| self.range_extractor_fwd.generate(kmers_fwd, rec_fwd, stats));

        stats.time_get_ranges += duration;
        let (duration, ranges_rev) =
            stage_time(|| self.range_extractor_rev.generate(kmers_rev, rec_rev, stats));
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
                let kmer_len = if range.vrange.has_header() {
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
                let kmer_len = if range.vrange.has_header() {
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
        // Taken out of the extractors (not copied) for the duration of the `&mut self` call below,
        // and put back afterwards so their capacity is reused: two `to_vec()` per pair were a
        // visible share of the allocator time.
        // --adaptive: a pair with a mate none of whose index ranges has at most
        // --adaptive-multiplicity positions -- every core of it a repeat of hundreds or thousands of
        // copies -- is the kind seeding misses the true copy of and the alignment window ranks it
        // out of: it gets more flank ties kept per core and more candidates aligned. The others
        // keep the defaults, so the effort stays with the reads that need it.
        let lim = self.options.args.adaptive_multiplicity;
        let hard = |ranges: &[crate::align::process::range_extractor::Range<F>]| !ranges.iter().any(|r| r.vrange.len() <= lim);
        let difficult = self.options.args.adaptive && (hard(ranges_fwd) || hard(ranges_rev));
        // The fewest copies of any core a mate kept (0: none kept).
        let minmult = |ranges: &[crate::align::process::range_extractor::Range<F>]| ranges.iter().map(|r| r.vrange.len()).min().unwrap_or(0);
        let (minmult_fwd, minmult_rev) = (minmult(ranges_fwd), minmult(ranges_rev));
        // --dump-candidates: each mate's ranges, their fewest copies, and those of at most 16 / 256.
        let range_stats = |ranges: &[crate::align::process::range_extractor::Range<F>]| {
            format!("{}\t{}\t{}\t{}", ranges.len(), minmult(ranges),
                ranges.iter().filter(|r| r.vrange.len() <= 16).count(), ranges.iter().filter(|r| r.vrange.len() <= 256).count())
        };
        let dumping = dump_this(rec_fwd.head());
        let ranges_dump = if dumping { Some((range_stats(ranges_fwd), range_stats(ranges_rev))) } else { None };
        // --mapq-copies: a pair whose mates' least repeated cores have many copies is placed among
        // copies seeding did not all reach. Its MAPQ ceiling falls with the log of both counts (a
        // mate with no core kept counts as --max-range-size copies). Human, cap 4096: a pair whose
        // best mate has a core of at most 16 copies is wrong 0.3-4 % of the time, one whose best
        // has more than 256 about 20 %.
        if self.options.args.mapq_copies > 0 {
            let full = self.options.args.max_range_size.max(1);
            let l = |x: usize| ((if x == 0 { full } else { x }) as f64).log2();
            let (lo, hi) = (l(minmult_fwd).min(l(minmult_rev)), l(minmult_fwd).max(l(minmult_rev)));
            let q = self.options.args.mapq_copies as f64 - self.options.args.mapq_copies_lo * lo - self.options.args.mapq_copies_hi * hi;
            COPIES_CAP.with(|u| u.set(q.clamp(0.0, 60.0) as u8));
        }
        let default_flex = self.seed_extractor_fwd.max_best_flex();
        if difficult {
            stats.adaptive_pairs += 1;
            self.seed_extractor_fwd.set_max_best_flex(self.options.args.adaptive_flex);
            self.seed_extractor_rev.set_max_best_flex(self.options.args.adaptive_flex);
            ALIGN_TOP_Y.with(|u| u.set(self.options.args.adaptive_top_y));
        }
        let (duration, _) = stage_time(|| self.seed_extractor_fwd.generate(ranges_fwd, stats));
        stats.time_range_header += duration;
        let (duration, _) = stage_time(|| self.seed_extractor_rev.generate(ranges_rev, stats));
        stats.time_range_header += duration;
        // A mate every core of which was dropped for more equally good flank copies than even the
        // recovery pass keeps: its copies are many, and the true one is among them only if all
        // are looked at. --adaptive-ties keeps up to that many; the partner places it among them.
        let ties = self.options.args.adaptive_ties;
        if self.options.args.adaptive && ties > 0 {
            for fwd in [true, false] {
                let (se, ranges) = if fwd { (&mut self.seed_extractor_fwd, ranges_fwd) } else { (&mut self.seed_extractor_rev, ranges_rev) };
                let (seeded, tied) = se.last_ties();
                if seeded == 0 && tied > 0 {
                    stats.adaptive_tied += 1;
                    se.set_max_best_flex(ties);
                    let (duration, _) = stage_time(|| se.generate(ranges, stats));
                    stats.time_range_header += duration;
                    se.set_max_best_flex(default_flex);
                }
            }
        }
        if difficult {
            self.seed_extractor_fwd.set_max_best_flex(default_flex);
            self.seed_extractor_rev.set_max_best_flex(default_flex);
        }
        let seeds_fwd = self.seed_extractor_fwd.take_seeds();
        let seeds_rev = self.seed_extractor_rev.take_seeds();
        if let Some(truth) = TRACE_TRUTH.get() {
            let flex = if difficult { self.options.args.adaptive_flex } else { self.options.args.max_best_flex };
            let lines = [
                trace_mate(self.db, truth, rec_fwd, "1", n_kmers_fwd, self.range_extractor_fwd.retrieve(), &seeds_fwd, flex, self.options.args.ranges),
                trace_mate(self.db, truth, rec_rev, "2", n_kmers_rev, self.range_extractor_rev.retrieve(), &seeds_rev, flex, self.options.args.ranges),
            ];
            if let Some((w, _)) = DUMP.get() {
                use std::io::Write as _;
                if let Some(w) = w.lock().unwrap().as_mut() {
                    for l in lines.into_iter().flatten() {
                        let _ = w.write_all(l.as_bytes());
                    }
                }
            }
        }
        // --mapq-copies-flank: the copy counts of the cap after the flank filter. A core with one
        // copy whose flanks match the read exactly makes the mate unique (human: 0.07 % of such
        // mates misplaced, against 2 % for a least repeated core of 2-16 copies); else the fewest
        // copies kept for a core.
        if self.options.args.mapq_copies > 0 && self.options.args.mapq_copies_flank {
            let full = self.options.args.max_range_size.max(1);
            let copies = |seeds: &[Seed]| {
                let (kept, exact) = flank_copies(seeds);
                if exact == 1 { 1 } else if kept == 0 { full } else { kept }
            };
            let l = |x: usize| (x as f64).log2();
            let (a, b) = (l(copies(&seeds_fwd)), l(copies(&seeds_rev)));
            let q = self.options.args.mapq_copies as f64 - self.options.args.mapq_copies_lo * a.min(b) - self.options.args.mapq_copies_hi * a.max(b);
            COPIES_CAP.with(|u| u.set(q.clamp(0.0, 60.0) as u8));
        }
        stats.seeds += seeds_fwd.len() + seeds_rev.len();

        // Seedless: no seed at all -- or, with --seedless-multiplicity N, none from a core of at
        // most N positions (with a raised --max-range-size, a mate seeded only by cores of
        // thousands of copies is placed among the copies its seeds reached).
        let lim = self.options.args.seedless_multiplicity;
        let seedless_mate = |seeds: &[Seed], re: &RE| seeds.is_empty() || (lim > 0 && re.ranges_at_most(lim) == 0);
        SEEDLESS.with(|u| u.set(seedless_mate(&seeds_fwd, &self.range_extractor_fwd) as u8 + seedless_mate(&seeds_rev, &self.range_extractor_rev) as u8));
        if std::env::var_os("FA_READ_STATS").is_some() {
            // Per pair: each mate's fewest core copies, and (ranges seeded, ranges dropped for
            // flank ties) of each mate.
            let (ta, tb) = (self.seed_extractor_fwd.last_ties(), self.seed_extractor_rev.last_ties());
            eprintln!("READSTAT\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                String::from_utf8_lossy(rec_fwd.head()).split_whitespace().next().unwrap_or("").trim_start_matches('@').split('/').next().unwrap_or(""),
                minmult_fwd, minmult_rev, ta.0, ta.1, tb.0, tb.1);
        }
        if let Some((rf, rr)) = ranges_dump {
            let (ta, tb) = (self.seed_extractor_fwd.last_ties(), self.seed_extractor_rev.last_ties());
            DUMP_READ.with(|d| *d.borrow_mut() = format!(
                "{rf}\t{}\t{}\t{}\t{}\t{rr}\t{}\t{}\t{}\t{}\t{}\t{}",
                seed_stats(&seeds_fwd), ta.0, ta.1, COPIES_CAP.with(|u| u.get()),
                seed_stats(&seeds_rev), tb.0, tb.1, difficult as u8, minmult_fwd.min(minmult_rev), minmult_fwd.max(minmult_rev)));
        }
        self.align_from_seeds(rec_fwd, rec_rev, &seeds_fwd, &seeds_rev, stats);
        DUMP_READ.with(|d| d.borrow_mut().clear());
        SEEDLESS.with(|u| u.set(0));
        ALIGN_TOP_Y.with(|u| u.set(0));
        COPIES_CAP.with(|u| u.set(60));

        self.seed_extractor_fwd.put_seeds(seeds_fwd);
        self.seed_extractor_rev.put_seeds(seeds_rev);
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
        if STAGE_TRACE {
            let started = std::time::Instant::now();
            stats.trace_counts.pairs_seeded += 1;
            self.align_from_seeds_traced(rec_fwd, rec_rev, seeds_fwd, seeds_rev, stats);
            stats.trace_counts.time_align_half += started.elapsed();
        } else {
            self.align_from_seeds_traced(rec_fwd, rec_rev, seeds_fwd, seeds_rev, stats);
        }
    }

    fn align_from_seeds_traced(
        &mut self,
        rec_fwd: &RefFastqRecord,
        rec_rev: &RefFastqRecord,
        seeds_fwd: &[Seed],
        seeds_rev: &[Seed],
        stats: &mut Stats,
    ) {
        // A read (pair) starts: no reference slice of the previous one is used again.
        crate::database::refseq2::reset_scratch();
        let dump_db = self.db;
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

        let (duration, anchors) = stage_time(|| {
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

        if STAGE_TRACE {
            let tc = &mut stats.trace_counts;
            tc.anchor_pairs += anchors.len();
            tc.pairs_no_anchor += anchors.is_empty() as usize;
            for AnchorPair(a1, a2) in anchors.iter() {
                tc.anchor_mates += a1.is_some() as usize + a2.is_some() as usize;
                tc.anchor_full_pairs += (a1.is_some() && a2.is_some()) as usize;
            }
        }

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

        let (duration, _) = stage_time(|| {
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
        //
        // TIES CROSS THAT BOUNDARY. The cut is a count, but what it cuts is a score order, and
        // equal keys keep arrival order -- which is index-lookup order, not evidence. A read whose
        // true locus ties the Zth anchor is therefore kept or dropped by list position alone. So
        // extend the window while the next anchor scores EXACTLY the same as the last one taken:
        // a tie is then settled by extension, which compares reference sequence, instead of by
        // where the index happened to report it. The window ends at the first strictly worse anchor
        // -- OR at `tie_cap`, whichever comes first. When the cap binds, the window stops mid-tie
        // and the remainder of that run is still cut by list order: the coin flip moves from Z to
        // the cap rather than disappearing.
        //
        // The key must match the extractor's sort key (`-(score(a1) + score(a2))`, missing side 0)
        // or the run this walks is not the run the sort produced.
        // AND THE TIE RUN IS CAPPED -- as a bound on a pathology that has NOT been observed, not as
        // a fix for a measured one. The concern is real in principle: with `-z` at 32 against ~15.7
        // anchors/read the loop is a no-op on ordinary reads and fires only where >32 anchors
        // already exist and score EQUALLY (a conserved marker hit by a family of near-identical
        // references), nothing bounds `anchors.len()` tightly there (headerless ranges emit up to
        // --max-range-size positions each without spending --ranges budget, and pairing is a
        // cross-product per reference), and the extender walks the whole slice doing a reference
        // fetch plus a full-read Hamming pass per entry. But measurement did not find it: on protal
        // an effectively-uncapped run (cap 64) was no slower than cap 2 and produced 4 more
        // successful alignments out of ~165,000. Keep the bound, do not claim it buys speed.
        let pair_anchor_score = |AnchorPair(a1, a2): &AnchorPair| {
            a1.as_ref().map_or(0, StdAnchorScore::score)
                + a2.as_ref().map_or(0, StdAnchorScore::score)
        };
        let screen_z = min(self.options.args.extend_top_z.max(1), anchors.len());
        let tie_cap = min(
            screen_z.saturating_mul(self.options.args.extend_tie_cap.max(1)),
            anchors.len(),
        );
        let mut extend_top_z = screen_z;
        let boundary_score = pair_anchor_score(&anchors[extend_top_z - 1]);
        while extend_top_z < tie_cap && pair_anchor_score(&anchors[extend_top_z]) == boundary_score {
            extend_top_z += 1;
        }
        debug_assert!(
            anchors.windows(2).all(|w| pair_anchor_score(&w[0]) >= pair_anchor_score(&w[1])),
            "anchors must reach the extension cut sorted by the extractor's key; the tie run walks \
             that order and silently truncates if some producer re-orders or leaves them unsorted"
        );

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
            // Deliberately `screen_z`, the CONFIGURED window, not the tie-extended one: this counter
            // is compared across runs, and a denominator that widens with a read's tie structure
            // would make two runs' funnels incomparable for reasons unrelated to what they screened.
            if holds_truth(&anchors[0..screen_z]) {
                stats.funnel_true_in_screen += 1;
            }
            // Both windows, because they answer different questions: `screen_z` is fixed-width and
            // therefore comparable across runs, while `extend_top_z` is what extension actually saw
            // and is therefore the one that attributes a loss to this stage.
            if holds_truth(&anchors[0..extend_top_z]) {
                stats.funnel_true_in_extended += 1;
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
            stats.funnel_true_in_screen_mates += holds_truth_mates(&anchors[0..screen_z]);
        }
        let pair_bonus = self.options.args.proper_pair_bonus;
        let mut time_extend = std::time::Duration::ZERO;
        let (duration, extended_count) = stage_time(|| {
            let extend_started = STAGE_TRACE.then(std::time::Instant::now);
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
            if let Some(t) = extend_started {
                time_extend = t.elapsed();
            }
            let bonus = pair_bonus;
            crate::align::common::sort_pairs_by_key(&mut anchors[0..extended_count], &mut self.fast.sort_scratch, |AnchorPair(a1, a2)| {
                let s = a1.as_ref().map_or(0, |a| a.score) + a2.as_ref().map_or(0, |a| a.score);
                -((s + if a1.is_some() && a2.is_some() { bonus } else { 0 }) as i64)
            });

            // assert!(extension_anchors.iter().all(|pair| pair.validate().is_ok()));
            extended_count
        });

        // Assumes sorted anchors !!
        let extension_anchors = &mut anchors[0..extended_count];
        let rescue_started = STAGE_TRACE.then(std::time::Instant::now);
        if STAGE_TRACE {
            stats.trace_counts.extend_window += extend_top_z;
            stats.trace_counts.extended += extended_count;
            stats.trace_counts.time_extend += time_extend;
        }

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
        if let Some(t) = rescue_started {
            stats.trace_counts.time_mate_rescue += t.elapsed();
        }
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


        let top_y = match ALIGN_TOP_Y.with(|u| u.get()) { 0 => self.options.args.align_top_y, y => y };
        // Assumes sorted anchors !!
        let anchors_len: usize = extension_anchors.len();
        let mut alignment_anchors =
            &mut extension_anchors[0..min(top_y, anchors_len)];
        if STAGE_TRACE {
            stats.trace_counts.align_window_pairs += alignment_anchors.len();
            for AnchorPair(a1, a2) in alignment_anchors.iter() {
                stats.trace_counts.align_window_mates += a1.is_some() as usize + a2.is_some() as usize;
            }
        }

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

        let (duration, _) = stage_time(|| {
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
            // Which mates the pair behind `best_pair_cost` has. A mate is only held to that cost when
            // the best pair already places it: when the best pair is one-sided, the other mate's
            // alignment in a later candidate is that mate's only chance to be written (see the lone
            // mate fallback below), so it keeps its full --min-ani budget.
            let mut best_pair_has = [false; 2];
            let unpaired_penalty = self.options.args.unpaired_penalty;
            let mut trace_total = AlignTrace::default();
            let mut stage_hist = [0usize; 6];
            let mut stage_hist_true = [0usize; 6];

            alignment_anchors
                .iter_mut()
                .enumerate()
                .for_each(|(_i, AnchorPair(a1, a2) )| {
                    let reference = match a1 {
                        Some(a) => self.db.get_reference(a.reference as usize).unwrap(),
                        None => self
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
                                let floor = supported_min_ani(a, min_ani, &self.options.args);
                                min_score_1 = Some(tighten_budget(
                                    ani_abort_score(floor, penalties::MISMATCH, alignable).abs(),
                                    // This pair can only win by costing less than the best pair; its
                                    // other mate costs at least 0, or the penalty when it is absent.
                                    best_pair_cost
                                        .filter(|_| best_pair_has[0])
                                        .map(|c| (c - if a2.is_none() { unpaired_penalty } else { 0 }).max(0)),
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
                                // let (duration, (score, cigar, status)) = stage_time(|| self.align.align(&query[qr], &reference[rr]));

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
                                let floor = supported_min_ani(a, min_ani, &self.options.args);
                                min_score_2 = Some(tighten_budget(
                                    ani_abort_score(floor, penalties::MISMATCH, alignable).abs(),
                                    // Mate 1 of this pair is already priced (ranked as the penalty
                                    // when absent or dropped), so what is left of the best pair's
                                    // cost is all mate 2 may spend.
                                    best_pair_cost
                                        .filter(|_| best_pair_has[1])
                                        .map(|c| (c - usable(a1).map_or(unpaired_penalty, |x| -x.score)).max(0)),
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
                                // let (duration, (score, cigar, status)) = stage_time(|| self.align.align(&query[qr], &reference[rr]));

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
                    let s1 = a1.as_ref().map_or(-unpaired_penalty, |a| a.score);
                    let s2 = a2.as_ref().map_or(-unpaired_penalty, |a| a.score);
                    if s1 > SUNK_SCORE && s2 > SUNK_SCORE {
                        let pair_cost = -(s1 + s2);
                        if pair_cost >= 0 && best_pair_cost.map_or(true, |c| pair_cost < c) {
                            best_pair_cost = Some(pair_cost);
                            best_pair_has = [a1.is_some(), a2.is_some()];
                        }
                    }
                });

            // --complete-singles: a candidate with one usable mate costs that mate plus the unpaired
            // penalty below -- about three mismatches -- so a wrong full pair of a few mismatches
            // beat a true single whose partner (unseeded: a repeat, an indel) lies right next to
            // it. Before ranking, rescue the partner next to each single that could still beat the
            // best full pair, and rank it complete.
            if self.options.args.complete_singles {
                let cost = |a: &Anchor| (-a.score).max(0) as i64;
                let full_best = alignment_anchors
                    .iter()
                    .filter_map(|p| Some(cost(usable(&p.0)?) + cost(usable(&p.1)?)))
                    .min();
                for p in alignment_anchors.iter_mut() {
                    let fwd_placed = match (usable(&p.0).is_some(), usable(&p.1).is_some()) {
                        (true, false) => true,
                        (false, true) => false,
                        _ => continue,
                    };
                    let placed = if fwd_placed { p.0.as_ref().unwrap() } else { p.1.as_ref().unwrap() };
                    if full_best.is_some_and(|b| cost(placed) >= b) || placed.reference_cigar_range.is_empty() {
                        continue;
                    }
                    let mate_seq = match (fwd_placed, placed.forward) {
                        (true, true) => self.rec_rev_revc.seq(),
                        (true, false) => rec_rev.seq(),
                        (false, true) => self.rec_fwd_revc.seq(),
                        (false, false) => rec_fwd.seq(),
                    };
                    let found = rescue_overhang(
                        self.db, &mut self.align, placed, mate_seq, min_ani,
                        self.options.args.min_query_coverage_mate, self.options.args.mate_rescue_margin as i32,
                    );
                    if let Some(r) = found {
                        stats.overhang_rescue_found += 1;
                        if fwd_placed { p.1 = Some(r) } else { p.0 = Some(r) }
                    }
                }
            }

            // A mate the aligner dropped is ranked like a mate with no anchor: both cost
            // --unpaired-penalty. Before, a missing mate cost 0 and a dropped one 1_000_000, so a
            // pair holding nothing but one dropped anchor outranked a pair with a good alignment and
            // a dropped partner, and the good alignment was never written. A pair with no usable
            // mate at all goes last.
            crate::align::common::sort_pairs_by_key(&mut alignment_anchors, &mut self.fast.sort_scratch, |AnchorPair(a1, a2)| {
                let cost = |u: Option<&Anchor>| u.map_or(unpaired_penalty, |a| -a.score) as i64;
                match (usable(a1), usable(a2)) {
                    (None, None) => i64::MAX,
                    (u1, u2) => cost(u1) + cost(u2),
                }
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
            pseudo_mapq = StdPairedAnchorMAPQ::mapq_from(extension_anchors, true, self.options.args.unpaired_penalty);
        } else if self.options.args.tie_mapq_cap > 0 && pseudo_mapq == 0 && !extension_anchors.is_empty() {
            // The anchors cannot tell the best pair from the runner-up; the alignments may.
            pseudo_mapq = StdPairedAnchorMAPQ::mapq_from(extension_anchors, true, self.options.args.unpaired_penalty)
                .min(self.options.args.tie_mapq_cap);
        }
        // The lowest cap applied below, for a MAPQ recomputed after the output's mate rescue.
        let mut mapq_ceiling: u8 = 60;
        // Mates the index gave no seed at all (every core dropped): placed by a rescue only.
        let seedless = SEEDLESS.with(|u| u.get());
        let seedless_cap = match seedless { 2 => self.options.args.seedless_pair_mapq, 1 => self.options.args.seedless_mate_mapq, _ => 60 };
        pseudo_mapq = pseudo_mapq.min(seedless_cap);
        mapq_ceiling = mapq_ceiling.min(seedless_cap);
        // Not part of the ceiling: the rescue's runner-up search below still runs, and its MAPQ
        // is held under this cap too.
        let copies_cap = COPIES_CAP.with(|u| u.get());
        pseudo_mapq = pseudo_mapq.min(copies_cap);

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

        let _reference = self.db.get_reference(*reference_id as usize).unwrap();

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
        //             Some(a) => self.db.get_reference(a.reference as usize).unwrap(),
        //             None => self.db.get_reference(a2.as_ref().unwrap().reference as usize).unwrap(),
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

        // LONE MATE FALLBACK. Only the winning pair is written, so a mate it has no usable alignment
        // for -- one-sided, or dropped by the aligner -- used to go unmapped even when that mate
        // aligned cleanly in another candidate. Mates whose own candidates are on different markers
        // (or whose partner has none) never share a pair, and against minibwa this was a third of
        // the reads flexalign missed. Such a mate now takes its own best alignment in the window; it
        // is written as a lone mate, so --lone-mate-mapq-penalty applies to it.
        let aligned_window = &extension_anchors[..min(top_y, extension_anchors.len())];
        let fallback = !self.options.args.no_lone_mate_fallback;
        // The winning pair's mates first; the fallback only for a mate the rescue next to its
        // partner (below) does not place either -- a mate's best alignment elsewhere, taken before
        // looking next to its partner, put it far from a correctly placed partner.
        let mut fwd_anchor = usable(&best_extended_anchor_pair.0);
        let mut rev_anchor = usable(&best_extended_anchor_pair.1);
        let mut rescued: Option<Anchor> = None;
        let mut fwd_out = resolve_mate(fwd_anchor, self.db, rec_fwd, &self.rec_fwd_revc);
        let mut rev_out = resolve_mate(rev_anchor, self.db, rec_rev, &self.rec_rev_revc);
        // A placed mate whose partner has no candidate at all: look for the partner next to it,
        // letting it hang off the end of the marker (`rescue_overhang`). Only behind a mate that
        // passes the output gates on its own -- the rescued one is then held to the pair-aware
        // floor (--min-query-coverage-mate), which is exactly the evidence this supplies.
        if !self.options.args.no_overhang_rescue && (fwd_anchor.is_none() != rev_anchor.is_none()) {
            let min_aligned_bases = self.options.args.min_query_coverage;
            let require_e2e = !self.options.args.allow_partial;
            let placed_ok = |m: &MateOut| {
                m.ani >= anchor_min_ani(m.anchor, min_ani, &self.options.args) && m.aligned_bases >= min_aligned_bases && (!require_e2e || m.end_to_end)
            };
            let (placed, mate_seq) = match (fwd_out.as_ref(), rev_out.as_ref()) {
                (Some(m), None) if placed_ok(m) => {
                    (m.anchor, if m.anchor.forward { self.rec_rev_revc.seq() } else { rec_rev.seq() })
                }
                (None, Some(m)) if placed_ok(m) => {
                    (m.anchor, if m.anchor.forward { self.rec_fwd_revc.seq() } else { rec_fwd.seq() })
                }
                _ => (fwd_anchor.or(rev_anchor).unwrap(), &[][..]),
            };
            if !mate_seq.is_empty() {
                stats.overhang_rescue_attempted += 1;
                rescued = rescue_overhang(
                    self.db,
                    &mut self.align,
                    placed,
                    mate_seq,
                    min_ani,
                    self.options.args.min_query_coverage_mate,
                    self.options.args.mate_rescue_margin as i32,
                );
                if let Some(a) = rescued.as_ref() {
                    stats.overhang_rescue_found += 1;
                    if fwd_out.is_none() {
                        fwd_out = resolve_mate(Some(a), self.db, rec_fwd, &self.rec_fwd_revc);
                    } else {
                        rev_out = resolve_mate(Some(a), self.db, rec_rev, &self.rec_rev_revc);
                    }
                }
            }
        }
        if fallback && rescued.is_none() {
            if fwd_anchor.is_none() {
                fwd_anchor = best_lone_mate(aligned_window, |p| &p.0);
                fwd_out = resolve_mate(fwd_anchor, self.db, rec_fwd, &self.rec_fwd_revc);
            }
            if rev_anchor.is_none() {
                rev_anchor = best_lone_mate(aligned_window, |p| &p.1);
                rev_out = resolve_mate(rev_anchor, self.db, rec_rev, &self.rec_rev_revc);
            }
        }

        // RESCUE-AWARE MAPQ. The MAPQ above was taken over the candidates as anchored: when the best
        // has one mate only, a runner-up with one mate elsewhere that aligns as well makes it 0 --
        // although the best's partner, found next to it just now, makes it a proper pair, and the
        // runner-up's partner may be nowhere near the runner-up. Rescue the partner of each such
        // runner-up the same way (up to three), and take the MAPQ from the pairs' costs: each mate's
        // alignment cost, the unpaired penalty for a partner not found.
        // Not for a pair already capped as a repeat: its candidates are known to be incomplete.
        if self.options.args.rescue_mapq && mapq_ceiling == 60 {
            if let (Some(_), Some(placed)) = (rescued.as_ref(), fwd_anchor.or(rev_anchor)) {
                let placed_fwd = fwd_anchor.is_some();
                let cost = |a: &Anchor| (-(a.score as i64)).max(0);
                let unpaired = self.options.args.unpaired_penalty.max(0) as i64;
                let mm_cost = penalties::MISMATCH as i64;
                let mate_seq_of = |a: &Anchor| match (placed_fwd, a.forward) {
                    (true, true) => self.rec_rev_revc.seq(),
                    (true, false) => rec_rev.seq(),
                    (false, true) => self.rec_fwd_revc.seq(),
                    (false, false) => rec_fwd.seq(),
                };
                // Both pairs' partners measured the same way: the gap-free scan's mismatches.
                let partner_cost = |a: &Anchor| {
                    partner_scan_mismatches(self.db, a, mate_seq_of(a), min_ani).map_or(unpaired, |mm| mm as i64 * mm_cost)
                };
                let best_cost = cost(placed) + partner_cost(placed);
                let mut alt_cost = i64::MAX;
                let mut tried = 0;
                for p in extension_anchors.iter() {
                    let (a, other) = if placed_fwd { (&p.0, &p.1) } else { (&p.1, &p.0) };
                    let Some(a) = a.as_ref() else { continue };
                    if a.score <= -1_000_000 || a.reference_cigar_range.is_empty() {
                        continue;
                    }
                    let near = a.reference == placed.reference
                        && a.reference_cigar_range.start.abs_diff(placed.reference_cigar_range.start) < MAX_INSERT as usize;
                    if near {
                        continue;
                    }
                    let c = match other.as_ref() {
                        Some(o) if o.score > -1_000_000 => cost(a) + cost(o),
                        _ => {
                            if tried == 8 {
                                continue;
                            }
                            tried += 1;
                            cost(a) + partner_cost(a)
                        }
                    };
                    alt_cost = alt_cost.min(c);
                }
                let per_mm = match self.options.args.mapq_per_mismatch { 0 => 30, v => v } as i64;
                let q = if alt_cost == i64::MAX {
                    crate::align::common::unique_mapq(self.options.args.mapq_unique as f64, best_cost as f64).round() as i64
                } else {
                    // Whole mismatches only: a margin below one is a tie.
                    (crate::align::common::per_mismatch_mapq(per_mm as f64, best_cost as f64)
                        * ((alt_cost - best_cost).max(0) / mm_cost) as f64).round().min(60.0) as i64
                };
                let q = (q as f64).min(crate::align::common::fit_cap(best_cost as f64)).round() as i64;
                pseudo_mapq = (q as u8).min(mapq_ceiling).min(copies_cap);
            }
        }

        if dump_this(rec_fwd.head()) {
            dump_candidates(dump_db, rec_fwd, extension_anchors, fwd_anchor, rev_anchor, pseudo_mapq);
        }

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
            m.ani >= anchor_min_ani(m.anchor, min_ani, &self.options.args)
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
                && m.ani >= anchor_min_ani(m.anchor, min_ani, &self.options.args)
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
        if STAGE_TRACE {
            stats.trace_counts.mates_resolved += fwd_out.is_some() as usize + rev_out.is_some() as usize;
            stats.trace_counts.mates_written += fwd_mate.is_some() as usize + rev_mate.is_some() as usize;
        }

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
                    self.options.args.lone_mate_mapq_penalty,
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
                    self.options.args.lone_mate_mapq_penalty,
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

    #[test]
    fn mismatches_capped_counts_like_a_byte_loop() {
        let q: Vec<u8> = (0..150u32).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        for shift in 0..20usize {
            let mut r = q.clone();
            for k in (0..r.len()).step_by(shift + 3) {
                r[k] = b'N';
            }
            for len in [0usize, 5, 8, 20, 37, 150] {
                let naive = q[..len].iter().zip(&r[..len]).filter(|(a, b)| a != b).count() as u32;
                assert_eq!(mismatches_capped(&q[..len], &r[..len], u32::MAX), naive);
                // Capped: exact up to the cap, and past it once the cap is exceeded.
                let capped = mismatches_capped(&q[..len], &r[..len], 3);
                assert!(if naive <= 3 { capped == naive } else { capped > 3 });
            }
        }
    }

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


/// Touches the cache line at the start of `bytes` without waiting for it.
#[inline(always)]
fn prefetch_line(bytes: &[u8]) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
        _mm_prefetch(bytes.as_ptr() as *const i8, _MM_HINT_T0);
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = bytes;
    }
}

/// --dump-candidates: one P line for the pair (read-level fields from `run_pair`, the final
/// MAPQ, where its mates were written) and one C line per aligned candidate pair, best first.
/// Columns: see DUMP_HEADER.
fn dump_candidates<DB: FlexalignDatabase>(db: &DB, rec_fwd: &RefFastqRecord, cands: &[AnchorPair], fa: Option<&Anchor>, ra: Option<&Anchor>, mapq: u8) {
    use std::fmt::Write as _;
    let Some((file, _)) = DUMP.get() else { return };
    let name = String::from_utf8_lossy(rec_fwd.head()).split_whitespace().next().unwrap_or("").trim_start_matches('@').to_string();
    let name = name.strip_suffix("/1").unwrap_or(&name).to_string();
    let rname = |a: &Anchor| db.get_rname(a.reference as usize).unwrap_or("?").to_string();
    let at = |a: Option<&Anchor>| a.filter(|a| !a.reference_cigar_range.is_empty())
        .map_or("-".to_string(), |a| format!("{}:{}", rname(a), a.reference_cigar_range.start));
    let mate = |a: &Option<Anchor>| -> String {
        match a.as_ref() {
            None => "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-".to_string(),
            // A candidate the window did not reach (or the aligner gave up on): its seeds only.
            Some(a) if a.score <= -1_000_000 || a.cigar.is_none() => format!(
                "{}\t{}\t-\t{}\t{}\t{}\t{}\t{}\t-\t-\t-\t-\t-\t-",
                rname(a), a.seeds.first().map_or(0, |s| s.rpos), if a.forward { '+' } else { '-' },
                a.seed_count, a.counter1, a.seeds.iter().map(|s| s.length).sum::<u32>(), a.mismatches),
            Some(a) => {
                let (mut m, mut x, mut i, mut d, mut s) = (0usize, 0usize, 0usize, 0usize, 0usize);
                if let Some(c) = a.cigar.as_ref() {
                    for &op in &c.0 {
                        match op { b'M' | b'=' => m += 1, b'X' => x += 1, b'I' => i += 1, b'D' => d += 1, b'S' => s += 1, _ => {} }
                    }
                }
                let ident = if m + x + i + d > 0 { m as f64 / (m + x + i + d) as f64 } else { 0.0 };
                format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
                    rname(a), a.reference_cigar_range.start, a.reference_cigar_range.end, if a.forward { '+' } else { '-' },
                    a.seed_count, a.counter1, a.seeds.iter().map(|s| s.length).sum::<u32>(), a.mismatches,
                    -(a.score as i64), m, x, i + d, s, ident)
            }
        }
    };
    let same = |a: &Option<Anchor>, b: Option<&Anchor>| match (a.as_ref(), b) {
        (Some(a), Some(b)) => a.reference == b.reference && a.reference_cigar_range == b.reference_cigar_range,
        (None, None) => true,
        _ => false,
    };
    let read = DUMP_READ.with(|d| d.borrow().clone());
    let mut out = String::new();
    let _ = writeln!(out, "P\t{name}\t{mapq}\t{}\t{}\t{}\t{}", cands.len(), at(fa), at(ra), if read.is_empty() { "-" } else { &read });
    for (k, p) in cands.iter().enumerate() {
        if p.0.is_none() && p.1.is_none() {
            continue;
        }
        let _ = writeln!(out, "C\t{name}\t{k}\t{}\t{}\t{}", (same(&p.0, fa) && same(&p.1, ra)) as u8, mate(&p.0), mate(&p.1));
    }
    use std::io::Write as _;
    if let Some(w) = file.lock().unwrap().as_mut() {
        let _ = w.write_all(out.as_bytes());
    }
}


/// The first line of a --dump-candidates file: the columns of its P and C lines.
pub const DUMP_HEADER: &str = "#P\tread\tmapq\tcandidates\twritten_fwd\twritten_rev\t\
f_ranges\tf_min_copies\tf_ranges_le16\tf_ranges_le256\tf_seeds\tf_exact_seeds\tf_cores_seeded\tf_min_kept\tf_min_exact\tf_ranges_seeded\tf_ranges_tie_dropped\tcopies_cap\t\
r_ranges\tr_min_copies\tr_ranges_le16\tr_ranges_le256\tr_seeds\tr_exact_seeds\tr_cores_seeded\tr_min_kept\tr_min_exact\tr_ranges_seeded\tr_ranges_tie_dropped\t\
adaptive\tmin_copies_lo\tmin_copies_hi\n\
#T\tread\tmate\tkmers\tranges\tranges_with_truth\tseeded\tflank_worse\ttie_cut\tbudget\tother\tseeds_at_truth\tmin_ties_cut\ttie_limit\tcut_flank_dist\tcut_truth_read_mm\tcut_best_read_mm\tcut_within0\tcut_within2\tcut_within5\n\
#F\tread\tmapq\twritten_fwd\twritten_rev\tf_mm\tr_mm\tf_cands\tr_cands\tvotes\trunner_up_votes\talt_pair_mm\talt_lone_mm\n\
#C\tread\trank\twritten\t\
f_ref\tf_start\tf_end\tf_strand\tf_seeds\tf_seeds_not_index\tf_seed_bases\tf_flank_mm\tf_cost\tf_match\tf_mismatch\tf_indel\tf_clip\tf_identity\t\
r_ref\tr_start\tr_end\tr_strand\tr_seeds\tr_seeds_not_index\tr_seed_bases\tr_flank_mm\tr_cost\tr_match\tr_mismatch\tr_indel\tr_clip\tr_identity\n";

/// `--dump-candidates`: the file every dumped pair is written to, and every how many pairs one is.
pub static DUMP: std::sync::OnceLock<(std::sync::Mutex<Option<Box<dyn std::io::Write + Send>>>, u64)> = std::sync::OnceLock::new();

/// Ends the --dump-candidates file (a .gz one gets its trailer) when dropped: held for the run.
pub struct DumpFinish;
impl Drop for DumpFinish {
    fn drop(&mut self) {
        if let Some((w, _)) = DUMP.get() {
            if let Some(mut w) = w.lock().unwrap().take() {
                let _ = w.flush();
            }
        }
    }
}

/// Whether the pair named `head` is dumped: one in --dump-every, by a hash of its name (the same
/// pairs whatever the thread count).
fn dump_this(head: &[u8]) -> bool {
    let Some((_, every)) = DUMP.get() else { return false };
    let name = head.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(head);
    let name = name.strip_suffix(b"/1").or_else(|| name.strip_suffix(b"/2")).unwrap_or(name);
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in name {
        h = (h ^ b as u64).wrapping_mul(0x100000001b3);
    }
    h % every.max(&1) == 0
}

/// `FA_TRACE_TRUTH=<file>` ("read\tmate\tcontig\tpos0" lines): the true placement of the reads to
/// trace; with --dump-candidates each traced mate gets a T line (`trace_mate`).
pub static TRACE_TRUTH: std::sync::OnceLock<std::collections::HashMap<String, (String, u64)>> = std::sync::OnceLock::new();

/// A T line: what became of the true locus at seeding. Of the mate's k-mers, how many found a range
/// (the rest are missing from the index: absent, or more copies than the cap); of the ranges, how
/// many contain the true copy, and of those how many seeded it, lost it to a copy with a better
/// flank, cut it at the flank-tie limit (`flex`), or were not walked (seed budget `max_ranges`
/// spent); and how many seeds landed at the true locus.
fn trace_mate<DB: FlexalignDatabase, const F: usize>(
    db: &DB, truth: &std::collections::HashMap<String, (String, u64)>, rec: &RefFastqRecord, mate: &str,
    n_kmers: usize, ranges: &[crate::align::process::range_extractor::Range<F>], seeds: &[Seed], flex: usize, max_ranges: usize,
) -> Option<String> {
    let name = String::from_utf8_lossy(rec.head()).split_whitespace().next().unwrap_or("").trim_start_matches('@').to_string();
    let name = name.strip_suffix("/1").or_else(|| name.strip_suffix("/2")).unwrap_or(&name).to_string();
    let (contig, pos0) = truth.get(&format!("{name}\t{mate}"))?;
    // The loaded index may not keep a name -> id table: build one from its id -> name table, once.
    static RID_OF: std::sync::OnceLock<std::collections::HashMap<String, u64>> = std::sync::OnceLock::new();
    let rid_of = RID_OF.get_or_init(|| (0..).map_while(|i| db.get_rname(i).map(|n| (n.to_string(), i as u64))).collect());
    let rid = *rid_of.get(contig)?;
    let (len, k) = (rec.seq().len() as i64, 31i64);
    let near = |rpos: u64, qpos: usize| {
        let (r, p, q) = (rpos as i64, *pos0 as i64, qpos as i64);
        (r - (p + q)).abs() <= 8 || (r - (p + len - q - k)).abs() <= 8
    };
    let (mut with_truth, mut seeded, mut flank_worse, mut tie_cut, mut budget, mut other) = (0, 0, 0, 0, 0, 0);
    let mut headered_before = 0usize;
    let mut buf = Vec::new();
    // Over the ranges cut at the tie limit that hold the true copy: the fewest best-flank ties.
    let mut min_ties = usize::MAX;
    // The first such range: its best flank distance, and the whole read against each tied copy
    // (gap-free mismatches, the better strand): the true copy's, the best tied copy's, and how many
    // tied copies come within 0 / 2 / 5 mismatches of that best.
    let mut cut_info: Option<(u32, usize, usize, usize, usize, usize)> = None;
    let fwd_seq = rec.seq().to_vec();
    let rc_seq: Vec<u8> = fwd_seq.iter().rev().map(|&b| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x }).collect();
    let read_mm = |rval: u64, rpos: u64, qpos: usize| -> usize {
        let Some(r) = db.get_reference(rval as usize) else { return usize::MAX };
        let mut best = usize::MAX;
        for (seq, start) in [(&fwd_seq, rpos as i64 - qpos as i64), (&rc_seq, rpos as i64 - (len - qpos as i64 - k))] {
            if start < 0 || start as usize + seq.len() > r.len() { continue }
            let w = &r[start as usize..start as usize + seq.len()];
            best = best.min(seq.iter().zip(w.iter()).filter(|(a, b)| a != b).count());
        }
        best
    };
    for r in ranges {
        let n = r.vrange.len();
        let idx = (0..n).find(|&j| { let (v, p) = r.vrange.value(j); v == rid && near(p, r.qpos) });
        if let Some(j) = idx {
            with_truth += 1;
            let got = seeds.iter().any(|s| s.rval == rid && (s.qpos as usize == r.qpos || s.qpos as usize == r.qpos + F / 2) && near(s.rpos - if s.mismatch == 0 { 0 } else { (F / 2) as u64 }, r.qpos));
            if got {
                seeded += 1;
            } else if r.vrange.has_header() {
                let min = crate::align::process::seed_extractor::flank_distances(&r.vrange, r.fmer.0 as u32, &mut buf);
                let d = buf.get(j).copied().unwrap_or(0) as u32;
                let at_min = buf.iter().filter(|&&x| x as u32 == min).count();
                if d <= min && at_min > flex && cut_info.is_none() {
                    let mms: Vec<usize> = buf.iter().enumerate().filter(|(_, &x)| x as u32 == min)
                        .map(|(i, _)| { let (v, p) = r.vrange.value(i); read_mm(v, p, r.qpos) }).collect();
                    let (tv, tp) = r.vrange.value(j);
                    let truth_mm = read_mm(tv, tp, r.qpos);
                    let best = mms.iter().copied().min().unwrap_or(usize::MAX);
                    let within = |e: usize| mms.iter().filter(|&&m| m <= best.saturating_add(e)).count();
                    cut_info = Some((min, truth_mm, best, within(0), within(2), within(5)));
                }
                if d > min { flank_worse += 1 } else if at_min > flex { tie_cut += 1; min_ties = min_ties.min(at_min) } else if headered_before >= max_ranges { budget += 1 } else { other += 1 }
            } else {
                other += 1;
            }
        }
        if r.vrange.has_header() {
            headered_before += 1;
        }
    }
    let at_locus = seeds.iter().filter(|s| s.rval == rid && ((s.rpos as i64) - *pos0 as i64).abs() <= len + 16).count();
    Some(format!("T\t{name}\t{mate}\t{n_kmers}\t{}\t{with_truth}\t{seeded}\t{flank_worse}\t{tie_cut}\t{budget}\t{other}\t{at_locus}\t{}\t{flex}\t{}\n", ranges.len(), if min_ties == usize::MAX { 0 } else { min_ties },
        cut_info.map_or("-\t-\t-\t-\t-\t-".to_string(), |c| format!("{}\t{}\t{}\t{}\t{}\t{}", c.0, c.1, c.2, c.3, c.4, c.5))))
}

/// The fewest seeds kept for a seeded core (the copies left after the flank filter) and the fewest
/// exact-flank seeds of a core with any (0: none), over a mate's seeds.
fn flank_copies(seeds: &[Seed]) -> (usize, usize) {
    // (k-mer start, exact): a core-only seed starts F/2 bases into its k-mer. Sorted, one run a core.
    let mut keys: Vec<(u32, bool)> = seeds.iter().map(|s| (s.qpos - if s.mismatch == 0 { 0 } else { 8 }, s.mismatch == 0)).collect();
    keys.sort_unstable_by_key(|k| k.0);
    let (mut min_kept, mut min_exact) = (usize::MAX, usize::MAX);
    for run in keys.chunk_by(|a, b| a.0 == b.0) {
        min_kept = min_kept.min(run.len());
        let exact = run.iter().filter(|k| k.1).count();
        if exact > 0 {
            min_exact = min_exact.min(exact);
        }
    }
    (if min_kept == usize::MAX { 0 } else { min_kept }, if min_exact == usize::MAX { 0 } else { min_exact })
}

/// `--min-ani-supported`: the identity floor of a candidate about to be aligned. One whose seeds cover
/// at least --min-ani-supported-bases read bases (more than one full k-mer) is held to the lower
/// floor and marked (`counter2`), so the output gate (`anchor_min_ani`) agrees with its budget.
/// Measured on other-strain reads against a marker DB: 94 % of the right placements below 0.95
/// identity had over 31 seed bases on their top candidate, against 1 % of the chance hits of reads
/// from outside every marker -- so the wider budget is spent on almost no off-target read.
fn supported_min_ani(a: &mut Anchor, min_ani: f64, args: &crate::options::Args) -> f64 {
    let low = args.min_ani_supported;
    if low > 0.0 && low < min_ani && a.seeds.iter().map(|s| s.length).sum::<u32>() >= args.min_ani_supported_bases {
        a.counter2 = 1;
        low
    } else {
        min_ani
    }
}

/// The identity floor an aligned anchor is held to at output: as `supported_min_ani` budgeted it.
fn anchor_min_ani(a: &Anchor, min_ani: f64, args: &crate::options::Args) -> f64 {
    if a.counter2 == 1 && args.min_ani_supported > 0.0 { args.min_ani_supported.min(min_ani) } else { min_ani }
}

/// A mate's seeds for --dump-candidates: count, those whose flanks match the read exactly, the
/// read positions (cores) seeded, and over those cores the fewest seeds kept (the copies left after
/// the flank filter) and the fewest exact-flank seeds of a core with any (0: none).
fn seed_stats(seeds: &[Seed]) -> String {
    let exact = seeds.iter().filter(|s| s.mismatch == 0).count();
    let mut per: std::collections::HashMap<u32, (usize, usize)> = std::collections::HashMap::new();
    for s in seeds {
        let e = per.entry(s.qpos - if s.mismatch == 0 { 0 } else { 8 }).or_default();
        e.0 += 1;
        e.1 += (s.mismatch == 0) as usize;
    }
    let min_kept = per.values().map(|v| v.0).min().unwrap_or(0);
    let min_exact = per.values().filter(|v| v.1 > 0).map(|v| v.1).min().unwrap_or(0);
    format!("{}\t{}\t{}\t{}\t{}", seeds.len(), exact, per.len(), min_kept, min_exact)
}

thread_local! {
    /// --dump-candidates: the read-level fields of the pair being aligned, set by `run_pair`.
    static DUMP_READ: std::cell::RefCell<String> = const { std::cell::RefCell::new(String::new()) };
    /// Mates of the pair being aligned without a seed from the index (0..=2).
    static SEEDLESS: std::cell::Cell<u8> = const { std::cell::Cell::new(0) };
    /// `--adaptive`: the alignment window of the pair being aligned, when not --align-top-y (0).
    static ALIGN_TOP_Y: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    /// `--mapq-copies`: the MAPQ ceiling of the pair being aligned from its cores' copy counts.
    static COPIES_CAP: std::cell::Cell<u8> = const { std::cell::Cell::new(60) };
}

