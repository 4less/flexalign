//! Microbenchmark: cost of flank-header resolution as a function of the number of elements N in a
//! range's header (ProjectShard.md flex filter / `Range::resolve_flex`).
//!
//! "range headers" in the pipeline is seed extraction, whose inner loop -- for each headered range
//! -- scans that range's `N` flank headers computing `HeaderSeq::dist` once each, finds the minimum,
//! then selects the cells within `--flank-slack` of it. This bench isolates that work and sweeps N,
//! so you can see whether the cost is linear in N and what one header element costs -- separately in
//! two regimes:
//!
//!   * hot      -- one range reused K times, so its `N*4` header bytes stay in L1/L2. This is the
//!                 pure ALU/cache cost (what you pay once the index is warm).
//!   * scattered-- many ranges laid out across a >LLC arena and visited in random order, so each
//!                 range's headers are a cache miss. This is the memory-latency regime -- the
//!                 in-RAM analogue of the major-page-fault cost that dominates on a cold multi-GB
//!                 index (there each miss is a *disk* read instead of a DRAM fetch).
//!
//! The resolution logic mirrors `Range::resolve_flex` exactly (two phases: cache each distance in a
//! `u8`, then select `dist <= min + slack`); `HeaderSeq::dist` is replicated inline (its u32 field
//! has no public constructor) -- it is a stable 5-op bit-parallel hamming distance.
//!
//! Run:  cargo +nightly run --release --example bench_header_dist [-- <flank_slack>]

use std::time::Instant;

/// Replica of `flexmap::values::HeaderSeq::dist`: hamming distance over F=16 2-bit flank bases.
#[inline(always)]
fn dist(header: u32, flex: u32) -> u32 {
    let x = header ^ flex;
    let a = x & 0x5555_5555;
    let b = (x & 0xAAAA_AAAA) >> 1;
    (a | b).count_ones()
}

/// One flex resolution over a range's `headers`/`positions`, identical in shape to
/// `Range::resolve_flex`: one pass caches every distance, a second selects the within-slack set and
/// reads those positions. Returns a checksum so nothing is optimized away.
#[inline(always)]
fn resolve(
    headers: &[u32],
    positions: &[u64],
    flex: u32,
    slack: u32,
    dist_buf: &mut Vec<u8>,
) -> u64 {
    dist_buf.clear();
    let mut min_dist = u32::MAX;
    for &h in headers {
        let d = dist(h, flex);
        dist_buf.push(d as u8);
        if d < min_dist {
            min_dist = d;
        }
    }
    let threshold = min_dist + slack;
    let mut chk = 0u64;
    for (i, &d) in dist_buf.iter().enumerate() {
        if d as u32 <= threshold {
            chk ^= positions[i]; // survivors read their VCell, as emit_seeds does
        }
    }
    chk
}

/// Cheap deterministic PRNG (SplitMix64) -- no external dep.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

fn make_range(rng: &mut Rng, n: usize) -> (Vec<u32>, Vec<u64>) {
    let headers: Vec<u32> = (0..n).map(|_| rng.next() as u32).collect();
    let positions: Vec<u64> = (0..n).map(|_| rng.next() & ((1 << 60) - 1)).collect();
    (headers, positions)
}

fn main() {
    let slack: u32 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(0);

    const NS: &[usize] = &[1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128, 256, 512, 1024, 2048, 4096];
    const HOT_TOUCHES: u64 = 200_000_000; // header elements touched per N in the hot sweep
    const ARENA_BYTES: usize = 512 << 20; // scattered arena >> LLC

    let mut rng = Rng(0xDEAD_BEEF);
    let mut dist_buf: Vec<u8> = Vec::new();
    let mut sink = 0u64;

    println!("flank_slack = {slack}");
    println!(
        "{:>6}  {:>12}  {:>12}   {:>12}  {:>12}",
        "N", "hot ns/rng", "hot ns/elem", "scat ns/rng", "scat ns/elem"
    );

    for &n in NS {
        // ---- hot: one range reused, cache-resident ----
        let (headers, positions) = make_range(&mut rng, n);
        let iters = (HOT_TOUCHES / n as u64).max(1);
        // warm up
        for _ in 0..1000 {
            sink ^= resolve(&headers, &positions, 0x1234_5678, slack, &mut dist_buf);
        }
        let t = Instant::now();
        for k in 0..iters {
            sink ^= resolve(&headers, &positions, k as u32, slack, &mut dist_buf);
        }
        let hot = t.elapsed().as_secs_f64() / iters as f64; // seconds per range
        let hot_ns_rng = hot * 1e9;
        let hot_ns_elem = hot_ns_rng / n as f64;

        // ---- scattered: many ranges across a >LLC arena, random visit order ----
        let bytes_per = n * (4 + 8);
        let m = (ARENA_BYTES / bytes_per).max(64);
        let ranges: Vec<(Vec<u32>, Vec<u64>)> = (0..m).map(|_| make_range(&mut rng, n)).collect();
        let order: Vec<usize> = {
            let mut o: Vec<usize> = (0..m).collect();
            // Fisher-Yates so access order defeats the prefetcher.
            for i in (1..m).rev() {
                let j = (rng.next() % (i as u64 + 1)) as usize;
                o.swap(i, j);
            }
            o
        };
        let t = Instant::now();
        for &idx in &order {
            let (h, p) = &ranges[idx];
            sink ^= resolve(h, p, idx as u32, slack, &mut dist_buf);
        }
        let scat = t.elapsed().as_secs_f64() / m as f64;
        let scat_ns_rng = scat * 1e9;
        let scat_ns_elem = scat_ns_rng / n as f64;

        println!(
            "{:>6}  {:>12.1}  {:>12.3}   {:>12.1}  {:>12.3}",
            n, hot_ns_rng, hot_ns_elem, scat_ns_rng, scat_ns_elem
        );
    }

    eprintln!("checksum {sink}");
}
