# flexmap change proposals (NOT implemented)

`../flexmap` is off-limits for this task (see `TASK.md`). Anything that would require editing it
gets written here instead: problem, exact call sites, proposed change, expected effect, cost of the
index invalidation it forces.

Rule of thumb for what lands here: the **query** side of k-mer selection lives in flexalign
(`src/align/process/kmer_extractor.rs`), the **build** side lives in flexmap (`src/build.rs`), and
they must select identically. So any change to *which cores get indexed* is a flexmap change even
when the edit looks local.

---

## P1 — syncmer selector rework (from `TODO.md`)

**Status:** proposal only. Blocked by the no-flexmap rule; also needs `kmerrs`.

`TODO.md` in the repo root documents three items against `ClosedSyncmer::is_minimizer` in the
`kmerr` crate — (1) hash the s-mer before comparing, to kill the `AAAA = 0x00` A-rich selection
bias; (2) raise `S` from 4 to 5–6 to cut the ~20 % tie rate; (3) the homopolymer case. The shipped
selector is an **open** syncmer at offset t=0, not a closed one (the `|| min == last` arm is
commented out), so nominal density is 1/12, not 2/12.

**Why it is a flexmap change:** the identical `ClosedSyncmer<C,S,L>` runs build-side in flexmap
`src/build.rs` (canonicalisation at ~lines 100-104 / 162-166 / 258-261 / 360-364). Changing the
rule in `kmerr` alone silently desynchronises build and query unless both are rebuilt together.

**Cost:** invalidates every existing index (bacteroides ~2 min, protal ~30 min for 32 GB), and needs
a `FLEXMAP_BLOB_VERSION` bump so the rebuild-on-mismatch guard catches stale blobs — exactly the
failure mode that produced the 0.03 %-aligned regression.

**Expected effect:** unmeasured. Would need its own benchmark run against all three datasets.

---

## P2 — _(add as encountered)_
