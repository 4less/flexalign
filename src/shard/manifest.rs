//! The run manifest -- what every pass must agree on (ProjectShard.md §6, §8).
//!
//! A shard pass and the rejoin communicate through `(chunk_id, idx_in_chunk)` read ids, which are
//! only meaningful if every pass decomposed the FASTQ into **exactly the same chunks**. Chunk
//! boundaries are a pure function of the input bytes and `buffer_size`, so those two things are
//! correctness inputs, not tuning knobs: a `buffer_size` that differs between a shard pass and the
//! rejoin silently re-pairs evidence with the wrong reads. Nothing downstream would catch it --
//! the ids still look valid, they just mean something else.
//!
//! Hence hurdle 6: `buffer_size` must not stay the `2usize.pow(24)` literal inlined at the
//! `process_fastq.rs` call sites. It lives here, and `verify_compatible` is the assertion.
//!
//! The parameters that shaped the evidence (`max_best_flex` in particular -- §4 makes the shard
//! apply the flex filter itself) are recorded for the same reason: evidence written at flex=16 and
//! replayed by a rejoin that believes it was written at flex=32 is not detectably wrong either.

use std::path::{Path, PathBuf};

use super::db::ShardRange;

pub const MANIFEST_VERSION: u32 = 1;

#[derive(Debug, thiserror::Error)]
pub enum ManifestError {
    #[error("manifest version {got}, expected {want}")]
    Version { got: u32, want: u32 },
    #[error(
        "buffer_size {got} does not match the manifest's {want} -- chunk boundaries would differ \
         and evidence would re-pair with the wrong reads"
    )]
    BufferSize { got: usize, want: usize },
    #[error("{field} is {got}, but the evidence was written with {want}")]
    Param { field: &'static str, got: String, want: String },
    #[error("input {path} is {got} bytes, but the evidence was written against {want} bytes")]
    InputLen { path: PathBuf, got: u64, want: u64 },
    #[error("input {path} content changed since the evidence was written")]
    InputContent { path: PathBuf },
    #[error("shards do not cover the key space: {0}")]
    Coverage(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Json(#[from] serde_json::Error),
}

/// Identifies an input FASTQ well enough to catch it being swapped or regenerated between passes.
///
/// The hash covers the length plus a bounded prefix, not the whole file: a shard run reads tens of
/// GB per pass and rehashing it N+1 times would cost more than the passes. That catches the
/// realistic failure -- a different or re-simulated file -- but not a deliberate collision, and not
/// an edit confined to the middle of a large file. Prefer keeping the manifest next to the run.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct InputFile {
    pub path: PathBuf,
    pub len: u64,
    pub prefix_hash: u64,
}

const HASH_PREFIX_BYTES: u64 = 1 << 20;

fn fnv1a(bytes: &[u8], mut h: u64) -> u64 {
    for b in bytes {
        h ^= *b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    h
}

impl InputFile {
    pub fn probe(path: impl AsRef<Path>) -> Result<Self, ManifestError> {
        use std::io::Read;
        let path = path.as_ref();
        let len = std::fs::metadata(path)?.len();
        let mut f = std::fs::File::open(path)?;
        let mut buf = Vec::new();
        f.by_ref().take(HASH_PREFIX_BYTES).read_to_end(&mut buf)?;
        Ok(Self {
            path: path.to_path_buf(),
            len,
            prefix_hash: fnv1a(&buf, 0xcbf2_9ce4_8422_2325),
        })
    }

    fn verify(&self) -> Result<(), ManifestError> {
        let now = Self::probe(&self.path)?;
        if now.len != self.len {
            return Err(ManifestError::InputLen {
                path: self.path.clone(),
                got: now.len,
                want: self.len,
            });
        }
        if now.prefix_hash != self.prefix_hash {
            return Err(ManifestError::InputContent { path: self.path.clone() });
        }
        Ok(())
    }
}

/// The k-mer geometry the index was built with. Evidence cells are raw `VCell`s decoded with
/// `Seed::from_flexmer::<K, C, F>`, so a mismatch here decodes the payload into wrong positions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Geometry {
    pub k: usize,
    pub c: usize,
    pub f: usize,
    pub s: usize,
}

/// The seed-extraction parameters the shard passes applied. `max_best_flex` is the one §4 spends
/// the exactness budget on: the shard applies it, so the rejoin cannot re-decide it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct SeedParams {
    pub max_best_flex: usize,
    pub max_ranges: usize,
    pub min_ranges: usize,
    pub max_range_size: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Manifest {
    pub version: u32,
    /// The determinism constant (§6). Chunk boundaries depend on this and the input bytes alone.
    pub buffer_size: usize,
    /// Number of chunks phase 0 produced. Sizes every evidence file's index (§8).
    pub n_chunks: u64,
    pub key_space: u64,
    pub shards: Vec<ShardRange>,
    pub geometry: Geometry,
    pub seed_params: SeedParams,
    pub inputs: Vec<InputFile>,
}

impl Manifest {
    pub fn new(
        buffer_size: usize,
        n_chunks: u64,
        key_space: u64,
        shards: Vec<ShardRange>,
        geometry: Geometry,
        seed_params: SeedParams,
        inputs: Vec<InputFile>,
    ) -> Self {
        Self {
            version: MANIFEST_VERSION,
            buffer_size,
            n_chunks,
            key_space,
            shards,
            geometry,
            seed_params,
            inputs,
        }
    }

    pub fn n_shards(&self) -> usize {
        self.shards.len()
    }

    /// Every key must be owned by exactly one shard, or the run silently drops (or double-counts)
    /// whatever falls in the gap.
    pub fn verify_coverage(&self) -> Result<(), ManifestError> {
        if self.shards.is_empty() {
            return Err(ManifestError::Coverage("no shards".into()));
        }
        if self.shards[0].lo != 0 {
            return Err(ManifestError::Coverage(format!(
                "first shard starts at {}, not 0",
                self.shards[0].lo
            )));
        }
        let last = self.shards.last().unwrap().hi;
        if last != self.key_space {
            return Err(ManifestError::Coverage(format!(
                "last shard ends at {}, key space is {}",
                last, self.key_space
            )));
        }
        for w in self.shards.windows(2) {
            if w[0].hi != w[1].lo {
                return Err(ManifestError::Coverage(format!(
                    "gap or overlap between {}..{} and {}..{}",
                    w[0].lo, w[0].hi, w[1].lo, w[1].hi
                )));
            }
        }
        Ok(())
    }

    /// Asserts that a pass about to run agrees with the evidence already on disk. This is the
    /// hurdle-6 check; call it before the rejoin reads a single record.
    pub fn verify_compatible(
        &self,
        buffer_size: usize,
        geometry: &Geometry,
        seed_params: &SeedParams,
    ) -> Result<(), ManifestError> {
        if self.version != MANIFEST_VERSION {
            return Err(ManifestError::Version { got: self.version, want: MANIFEST_VERSION });
        }
        if self.buffer_size != buffer_size {
            return Err(ManifestError::BufferSize { got: buffer_size, want: self.buffer_size });
        }
        if &self.geometry != geometry {
            return Err(ManifestError::Param {
                field: "geometry",
                got: format!("{:?}", geometry),
                want: format!("{:?}", self.geometry),
            });
        }
        if &self.seed_params != seed_params {
            return Err(ManifestError::Param {
                field: "seed_params",
                got: format!("{:?}", seed_params),
                want: format!("{:?}", self.seed_params),
            });
        }
        self.verify_coverage()
    }

    /// Re-probes every input. Separate from `verify_compatible` because it touches the filesystem.
    pub fn verify_inputs(&self) -> Result<(), ManifestError> {
        self.inputs.iter().try_for_each(InputFile::verify)
    }

    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), ManifestError> {
        let json = serde_json::to_string_pretty(self)?;
        std::fs::write(path, json)?;
        Ok(())
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self, ManifestError> {
        let json = std::fs::read_to_string(path)?;
        let m: Manifest = serde_json::from_str(&json)?;
        if m.version != MANIFEST_VERSION {
            return Err(ManifestError::Version { got: m.version, want: MANIFEST_VERSION });
        }
        Ok(m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn geometry() -> Geometry {
        Geometry { k: 31, c: 15, f: 16, s: 4 }
    }

    fn params() -> SeedParams {
        SeedParams { max_best_flex: 16, max_ranges: 15, min_ranges: 4, max_range_size: 256 }
    }

    fn manifest(shards: Vec<ShardRange>, key_space: u64) -> Manifest {
        Manifest::new(1 << 24, 1800, key_space, shards, geometry(), params(), vec![])
    }

    fn tmp(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("flexalign-manifest-{}-{}", std::process::id(), name));
        p
    }

    #[test]
    fn roundtrips_through_json() {
        let path = tmp("roundtrip");
        let m = manifest(ShardRange::split_even(1 << 20, 10), 1 << 20);
        m.save(&path).unwrap();
        assert_eq!(Manifest::load(&path).unwrap(), m);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn a_mismatched_buffer_size_is_rejected() {
        // The hurdle-6 failure: chunk boundaries would differ, so read ids would mean different
        // reads in the rejoin than they did in the shard pass.
        let m = manifest(ShardRange::split_even(1 << 20, 4), 1 << 20);
        let err = m.verify_compatible(1 << 26, &geometry(), &params());
        assert!(matches!(err, Err(ManifestError::BufferSize { got, want })
            if got == 1 << 26 && want == 1 << 24));

        m.verify_compatible(1 << 24, &geometry(), &params()).expect("matching buffer_size");
    }

    #[test]
    fn a_retuned_flex_threshold_is_rejected() {
        // §4: the shard applies the flex filter, so evidence written at 16 is not evidence at 32.
        let m = manifest(ShardRange::split_even(1 << 20, 4), 1 << 20);
        let mut p = params();
        p.max_best_flex = 32;
        assert!(matches!(
            m.verify_compatible(1 << 24, &geometry(), &p),
            Err(ManifestError::Param { field: "seed_params", .. })
        ));
    }

    #[test]
    fn a_different_geometry_is_rejected() {
        let m = manifest(ShardRange::split_even(1 << 20, 4), 1 << 20);
        let mut g = geometry();
        g.f = 8;
        assert!(matches!(
            m.verify_compatible(1 << 24, &g, &params()),
            Err(ManifestError::Param { field: "geometry", .. })
        ));
    }

    #[test]
    fn split_even_shards_always_cover() {
        for n in [1usize, 2, 3, 10, 64] {
            let m = manifest(ShardRange::split_even(1 << 20, n), 1 << 20);
            m.verify_coverage().unwrap_or_else(|e| panic!("n={}: {}", n, e));
        }
    }

    #[test]
    fn gaps_and_short_coverage_are_rejected() {
        let hole = vec![ShardRange::new(0, 16), ShardRange::new(32, 64)];
        assert!(matches!(manifest(hole, 64).verify_coverage(), Err(ManifestError::Coverage(_))));

        let short = vec![ShardRange::new(0, 16)];
        assert!(matches!(manifest(short, 64).verify_coverage(), Err(ManifestError::Coverage(_))));

        let late = vec![ShardRange::new(16, 64)];
        assert!(matches!(manifest(late, 64).verify_coverage(), Err(ManifestError::Coverage(_))));

        assert!(matches!(manifest(vec![], 64).verify_coverage(), Err(ManifestError::Coverage(_))));
    }

    #[test]
    fn input_probe_detects_a_changed_file() {
        let path = tmp("input");
        std::fs::write(&path, b"@r1\nACGT\n+\nIIII\n").unwrap();
        let probe = InputFile::probe(&path).unwrap();
        probe.verify().expect("unchanged file verifies");

        // Same length, different bytes -- the length check alone would miss this.
        std::fs::write(&path, b"@r1\nTGCA\n+\nIIII\n").unwrap();
        assert!(matches!(probe.verify(), Err(ManifestError::InputContent { .. })));

        std::fs::write(&path, b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n").unwrap();
        assert!(matches!(probe.verify(), Err(ManifestError::InputLen { .. })));

        std::fs::remove_file(&path).ok();
    }
}
