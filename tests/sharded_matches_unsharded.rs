//! The sharded pipeline must give the unsharded run's alignments, record for record.
//!
//! The two are separate code paths around shared stages (seed walk, alignment half), and twice
//! they drifted apart unnoticed: the join built its alignment half without the anchor gate, and the
//! shard passes could not do the min-ranges recovery. This runs the real binary on a generated
//! marker set -- with a segment copied, mutated, into many markers, so ranges are repetitive enough
//! to be cut by --max-best-flex and to trigger the recovery -- unsharded and on 2 and 4 shards, for
//! both packed index formats, and compares the SAM records.
//!
//! Fast-hit is off on both sides: the sharded pipeline has no fast-hit path (a pair pinned by its
//! unique k-mers skips seeding and is placed with its own MAPQ), so where the unsharded run's probe
//! turns fast-hit on -- as it does here, every read coming from a marker -- the MAPQ of those pairs
//! differs. On marker-DB metagenomes the probe turns it off (few reads hit).

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4) as usize]
    }
}

fn mutate(rng: &mut Rng, s: &mut [u8], per_mille: u64) {
    for b in s.iter_mut() {
        if rng.below(1000) < per_mille {
            *b = rng.base();
        }
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// 300 markers of 0.8-2.5 kb named like protal's (`<tiid>_<gene>`); a 600 bp segment sits,
/// 2% mutated per copy, in 80 of them. 4000 read pairs of 150 bp, inserts of 300-400, 1% errors,
/// a quarter of them from the copied segment.
fn write_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let family: Vec<u8> = (0..600).map(|_| rng.base()).collect();
    let mut markers: Vec<Vec<u8>> = Vec::new();
    let mut family_at: Vec<Option<usize>> = Vec::new();
    for i in 0..300 {
        let len = 800 + rng.below(1700) as usize;
        let mut m: Vec<u8> = (0..len).map(|_| rng.base()).collect();
        if i % 4 == 0 && m.len() > 700 {
            let mut copy = family.clone();
            mutate(&mut rng, &mut copy, 20);
            let at = rng.below((m.len() - 600) as u64) as usize;
            m[at..at + 600].copy_from_slice(&copy);
            family_at.push(Some(at));
        } else {
            family_at.push(None);
        }
        markers.push(m);
    }
    let reference = dir.join("ref.fna");
    let mut f = fs::File::create(&reference).unwrap();
    for (i, m) in markers.iter().enumerate() {
        writeln!(f, ">{}_{}", 1000 + i / 3, i).unwrap();
        f.write_all(m).unwrap();
        writeln!(f).unwrap();
    }
    let (r1, r2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    let (mut f1, mut f2) = (fs::File::create(&r1).unwrap(), fs::File::create(&r2).unwrap());
    let qual = vec![b'I'; 150];
    let mut n = 0;
    while n < 4000 {
        let i = rng.below(markers.len() as u64) as usize;
        let m = &markers[i];
        let insert = 300 + rng.below(100) as usize;
        if m.len() < insert {
            continue;
        }
        let start = match family_at[i] {
            Some(at) if n % 4 == 0 => at.saturating_sub(100).min(m.len() - insert),
            _ => rng.below((m.len() - insert + 1) as u64) as usize,
        };
        let frag = &m[start..start + insert];
        let mut a = frag[..150].to_vec();
        let mut b = revcomp(&frag[insert - 150..]);
        mutate(&mut rng, &mut a, 10);
        mutate(&mut rng, &mut b, 10);
        for (f, s) in [(&mut f1, &a), (&mut f2, &b)] {
            writeln!(f, "@r{n}").unwrap();
            f.write_all(s).unwrap();
            writeln!(f, "\n+").unwrap();
            f.write_all(&qual).unwrap();
            writeln!(f).unwrap();
        }
        n += 1;
    }
    (reference, r1, r2)
}

fn run(reference: &Path, r1: &Path, r2: &Path, out: &Path, extra: &[&str]) -> Vec<String> {
    let status = Command::new(env!("CARGO_BIN_EXE_flexalign"))
        .arg("-r").arg(reference)
        .arg("-1").arg(r1)
        .arg("-2").arg(r2)
        .arg("-0").arg(out)
        .args(["-t", "4"])
        .args(extra)
        .output()
        .expect("run flexalign");
    assert!(status.status.success(), "flexalign {extra:?} failed:\n{}", String::from_utf8_lossy(&status.stderr));
    // The sharded run writes <out>.sam, the unsharded one <out>.
    let sam = if out.exists() { out.to_path_buf() } else { out.with_extension("sam") };
    let text = fs::read_to_string(&sam).unwrap_or_else(|e| panic!("{}: {e}", sam.display()));
    let mut records: Vec<String> = text.lines().filter(|l| !l.starts_with('@')).map(String::from).collect();
    records.sort();
    records
}

#[test]
fn sharded_runs_give_the_unsharded_alignments() {
    let dir = std::env::temp_dir().join(format!("flexalign-shard-parity-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let (reference, r1, r2) = write_fixture(&dir);
    for format in ["hashed"] {
        let base = run(&reference, &r1, &r2, &dir.join(format!("{format}-none")), &["--no-shards", "--no-fast-hit"]);
        assert!(base.len() > 4000, "{format}: only {} records unsharded", base.len());
        for n in ["2", "4"] {
            let sharded = run(&reference, &r1, &r2, &dir.join(format!("{format}-s{n}")), &["--shards", n, "--no-fast-hit"]);
            let differ = base.iter().zip(&sharded).filter(|(a, b)| a != b).count() + base.len().abs_diff(sharded.len());
            assert!(
                differ == 0,
                "{format}, {n} shards: {differ} of {} records differ from the unsharded run; first: {:?}",
                base.len(),
                base.iter().zip(&sharded).find(|(a, b)| a != b)
            );
        }
    }
    fs::remove_dir_all(&dir).unwrap();
}
