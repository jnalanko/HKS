//! Hierarchy-aware smoothing for HKS lookup output.
//!
//! Implements Algorithm S1 from the HKS paper: a two-phase scan that identifies
//! windows exhibiting a specific → general → specific pattern in the category
//! hierarchy and reassigns interior intervals to the LCA of the flanking anchors.

use std::collections::HashSet;
use std::io::{BufWriter, Read, Write};

use crate::lca_tree::LcaTree;

// ---------------------------------------------------------------------------
// Interval representation
// ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
    pub feature: usize,     // node ID in the hierarchy
    pub originally_none: bool, // true if this interval was "none" / unmatched in the input
}

// ---------------------------------------------------------------------------
// Core smoothing algorithm (port of Algorithm S1)
// ---------------------------------------------------------------------------

/// Smooth a single query's intervals in-place using the hierarchy.
/// Returns the number of feature reassignments made.
pub fn smooth_intervals(intervals: &mut Vec<Interval>, tree: &LcaTree, max_gap: u64) -> u64 {
    if intervals.len() < 2 {
        return 0;
    }
    let mut total_reassignments = 0u64;
    let n = intervals.len();

    loop {
        let mut changed = false;
        let mut was_related = vec![false; n];
        let mut i = 0;

        while i < n {
            let window_start = i;
            let mut related: Vec<usize> = vec![i];
            let mut first_unrelated: Option<usize> = None;

            // ------------------------------------------------------------------
            // Ascending phase: scan rightward accepting features that are
            // ancestors of last_rel_feat (i.e. more general / closer to root).
            // Unrelated features on other branches are skipped but their ancestor
            // paths are added to `disallowed` so we stop if we'd cross into them.
            // ------------------------------------------------------------------
            let mut disallowed: HashSet<usize> = HashSet::new();
            let mut last_rel_idx = i;
            let mut last_rel_feat = intervals[i].feature;
            let mut last_rel_end = intervals[i].end;
            let mut j = i;

            while j + 1 < n {
                let nf = intervals[j + 1].feature;
                let ns = intervals[j + 1].start;

                if ns > last_rel_end + max_gap {
                    break;
                }

                if tree.is_ancestor(nf, last_rel_feat) {
                    // nf IS an ancestor of last_rel_feat → nf is more general → extend
                    if disallowed.contains(&nf) {
                        break;
                    }
                    last_rel_feat = nf;
                    last_rel_end = intervals[j + 1].end;
                    last_rel_idx = j + 1;
                    related.push(j + 1);
                    j += 1;
                } else if !tree.is_ancestor(last_rel_feat, nf) {
                    // Neither is ancestor of the other → different branches
                    disallowed.extend(tree.ancestors(nf));
                    if first_unrelated.is_none() && !was_related[j + 1] {
                        first_unrelated = Some(j + 1);
                    }
                    j += 1;
                } else {
                    // last_rel_feat IS ancestor of nf → nf is more specific → stop
                    break;
                }
            }

            let mut window_end = last_rel_idx;

            // ------------------------------------------------------------------
            // Descending phase: continue rightward from the peak, accepting
            // features that are descendants of last_rel_feat (more specific).
            // ------------------------------------------------------------------
            let peak_feat = intervals[window_end].feature;
            let mut k = window_end;

            while k + 1 < n {
                let nf = intervals[k + 1].feature;
                let ns = intervals[k + 1].start;

                if ns > last_rel_end + max_gap {
                    break;
                }

                if tree.is_ancestor(last_rel_feat, nf) {
                    // last_rel_feat IS ancestor of nf → nf is more specific → extend
                    last_rel_feat = nf;
                    last_rel_end = intervals[k + 1].end;
                    last_rel_idx = k + 1;
                    related.push(k + 1);
                    k += 1;
                } else if !tree.is_ancestor(nf, last_rel_feat) {
                    // Neither is ancestor → unrelated
                    if tree.is_ancestor(peak_feat, nf) {
                        // nf is a descendant of peak → would restart a new ascending window → stop
                        break;
                    }
                    if first_unrelated.is_none() && !was_related[k + 1] {
                        first_unrelated = Some(k + 1);
                    }
                    k += 1;
                } else {
                    // nf IS ancestor of last_rel_feat → nf is more general → stop
                    break;
                }
            }

            window_end = last_rel_idx;

            // Drop the last element from related (boundary stays unchanged)
            if related.len() >= 2 {
                related.pop();
            }
            for &idx in &related {
                was_related[idx] = true;
            }

            // ------------------------------------------------------------------
            // Reassignment: replace interior features with LCA(left, right) when
            // they are strictly more general (shallower) than the LCA.
            // ------------------------------------------------------------------
            if window_end > window_start {
                let left = intervals[window_start].feature;
                let right = intervals[window_end].feature;
                let lca = tree.lca(left, right);
                for w in (window_start + 1)..window_end {
                    let orig = intervals[w].feature;
                    // is_ancestor(orig, lca) means orig IS an ancestor of lca,
                    // i.e. orig is more general than lca → replace with lca
                    if tree.is_ancestor(orig, lca) && orig != lca {
                        intervals[w].feature = lca;
                        intervals[w].originally_none = false;
                        changed = true;
                        total_reassignments += 1;
                    }
                }
            }

            // Advance i
            i = match first_unrelated {
                Some(fu) if fu > window_end && window_end > window_start => window_end,
                Some(fu) => fu,
                None if window_end > i => window_end,
                None => {
                    let mut next = i;
                    while next < n && was_related[next] {
                        next += 1;
                    }
                    next
                }
            };
        }

        if !changed {
            break;
        }
    }
    total_reassignments
}

// ---------------------------------------------------------------------------
// Merge adjacent contiguous intervals with the same feature
// ---------------------------------------------------------------------------

/// Merges adjacent intervals that have the same feature and are contiguous.
/// Returns (merged intervals, number of intervals eliminated).
pub fn merge_intervals(intervals: Vec<Interval>) -> (Vec<Interval>, u64) {
    if intervals.is_empty() {
        return (intervals, 0);
    }
    let n_in = intervals.len();
    let mut out: Vec<Interval> = Vec::with_capacity(n_in);
    let mut cur = intervals.into_iter();
    let mut current = cur.next().unwrap();
    for next in cur {
        if next.feature == current.feature
            && next.start == current.end
            && next.originally_none == current.originally_none
        {
            current.end = next.end;
        } else {
            out.push(current);
            current = next;
        }
    }
    out.push(current);
    let eliminated = (n_in - out.len()) as u64;
    (out, eliminated)
}

// ---------------------------------------------------------------------------
// Streaming smooth processor: parse TSV → smooth per query → write TSV
// ---------------------------------------------------------------------------

#[derive(Default)]
pub struct SmoothStats {
    pub reads_processed: u64,
    pub intervals_in: u64,
    pub intervals_smoothed: u64,
    pub intervals_merged: u64,
    pub intervals_out: u64,
}

/// Resolve a feature token from the TSV input into a node ID.
/// Returns (node_id, originally_none).
fn resolve_feature(
    token: &str,
    name_to_id: &std::collections::HashMap<String, usize>,
    root_id: usize,
    uses_names: bool,
) -> (usize, bool) {
    if uses_names {
        if token == "none" {
            (root_id, true)
        } else {
            let id = name_to_id.get(token)
                .unwrap_or_else(|| panic!("Unknown feature name in input: {:?}", token));
            (*id, false)
        }
    } else {
        // Numeric ID mode
        if token == "-" {
            (root_id, true)
        } else {
            let id: usize = token.parse()
                .unwrap_or_else(|_| panic!("Cannot parse feature ID: {:?}", token));
            (id, false)
        }
    }
}

/// Format a feature for output.
fn format_feature(
    feature: usize,
    originally_none: bool,
    names: &[String],
    root_id: usize,
    uses_names: bool,
) -> String {
    if originally_none && feature == root_id {
        // Was none and smoothing didn't resolve it → keep as none
        if uses_names { "none".to_string() } else { "-".to_string() }
    } else if uses_names {
        names[feature].to_string()
    } else {
        feature.to_string()
    }
}

impl SmoothStats {
    /// Fold one group's partial stats into this accumulator.
    fn merge(&mut self, o: &SmoothStats) {
        self.reads_processed += o.reads_processed;
        self.intervals_in += o.intervals_in;
        self.intervals_smoothed += o.intervals_smoothed;
        self.intervals_merged += o.intervals_merged;
        self.intervals_out += o.intervals_out;
    }
}

/// Smooth one query group given as a raw byte slice of the input TSV.
///
/// Every line in `slice` shares one query id. Parses the slice, runs
/// smooth + merge, and renders the result to a freshly allocated byte buffer.
/// Doing the parse *and* the output formatting here (not just the smoothing
/// compute) is what lets the parallel path actually scale: parsing and
/// formatting tens of millions of lines dominate the wall time, so both must
/// run on the worker thread. Parsing borrows byte subslices of the input —
/// no per-line `String` allocation, unlike `reader.lines()`.
fn smooth_group(
    slice: &[u8],
    tree: &LcaTree,
    names: &[String],
    name_to_id: &std::collections::HashMap<String, usize>,
    root_id: usize,
    uses_names: bool,
    max_gap: u64,
) -> (SmoothStats, Vec<u8>) {
    let mut query_id: &str = "";
    let mut have_qid = false;
    let mut intervals: Vec<Interval> = Vec::new();

    for raw in slice.split(|&b| b == b'\n') {
        let line = raw.trim_ascii();
        if line.is_empty() {
            continue;
        }
        let mut cols = line.splitn(4, |&b| b == b'\t');
        let qid_b = cols.next().expect("missing query column");
        let start_b = cols.next().expect("missing start coordinate");
        let end_b = cols.next().expect("missing end coordinate");
        let feat_b = cols.next().expect("missing feature column");

        if !have_qid {
            query_id = std::str::from_utf8(qid_b).expect("non-UTF8 query name");
            have_qid = true;
        }
        let start: u64 = std::str::from_utf8(start_b).ok()
            .and_then(|s| s.parse().ok())
            .expect("bad start coordinate");
        let end: u64 = std::str::from_utf8(end_b).ok()
            .and_then(|s| s.parse().ok())
            .expect("bad end coordinate");
        let feat_token = std::str::from_utf8(feat_b).expect("non-UTF8 feature name");

        let (feature, originally_none) = resolve_feature(feat_token, name_to_id, root_id, uses_names);
        intervals.push(Interval { start, end, feature, originally_none });
    }

    let n_in = intervals.len() as u64;
    let reassigned = smooth_intervals(&mut intervals, tree, max_gap);
    let (merged, eliminated) = merge_intervals(intervals);
    let n_out = merged.len() as u64;

    // Grow from empty rather than reserving slice.len(): smoothing+merging
    // often collapses a huge input group (e.g. a whole chromosome) to a handful
    // of intervals, so reserving the input size would mmap hundreds of MB per
    // worker for a few bytes of output — concurrent large reservations serialize
    // on the kernel's mmap lock and cause a hard regression at high thread counts.
    let mut out: Vec<u8> = Vec::new();
    for iv in &merged {
        let feat_str = format_feature(iv.feature, iv.originally_none, names, root_id, uses_names);
        writeln!(out, "{}\t{}\t{}\t{}", query_id, iv.start, iv.end, feat_str)
            .expect("formatting error");
    }

    let stats = SmoothStats {
        reads_processed: 1,
        intervals_in: n_in,
        intervals_smoothed: reassigned,
        intervals_merged: eliminated,
        intervals_out: n_out,
    };
    (stats, out)
}

// ---------------------------------------------------------------------------
// Streaming grouper: split the input byte stream into per-query groups without
// ever holding the whole file in memory.
// ---------------------------------------------------------------------------

/// Immutable configuration threaded through the streaming parallel path.
struct SmoothCfg<'a> {
    tree: &'a LcaTree,
    names: &'a [String],
    name_to_id: &'a std::collections::HashMap<String, usize>,
    root_id: usize,
    max_gap: u64,
    max_groups: usize, // dispatch a batch once this many groups have closed …
    max_bytes: usize,  // … or once the closed groups reach this many bytes.
}

/// Accumulates data lines into per-query byte groups.
///
/// `buf` lays out `[ closed groups … | currently-open group ]`. The open group
/// is the maximal run of trailing lines that still share `open_qid`; it cannot
/// be dispatched yet because a later block may extend it. Closed groups (all of
/// `buf[..open_start]`) are handed off in batches and their bytes reclaimed, so
/// peak memory is bounded by one batch plus the single largest query group —
/// never the whole input.
#[derive(Default)]
struct Grouper {
    buf: Vec<u8>,
    closed: Vec<(usize, usize)>, // ranges into `buf`, all within [0, open_start)
    open_start: usize,           // start of the open group; also == closed-byte count
    open_qid: Vec<u8>,
    have_open: bool,
}

impl Grouper {
    /// Append one trimmed, non-empty data line (no trailing newline), closing
    /// the previous group first if this line begins a new query id.
    fn push_line(&mut self, line: &[u8]) {
        let qid = match line.iter().position(|&b| b == b'\t') {
            Some(k) => &line[..k],
            None => line,
        };
        let ls = self.buf.len();
        if !self.have_open {
            self.open_start = ls;
            self.open_qid.clear();
            self.open_qid.extend_from_slice(qid);
            self.have_open = true;
        } else if self.open_qid.as_slice() != qid {
            self.closed.push((self.open_start, ls));
            self.open_start = ls;
            self.open_qid.clear();
            self.open_qid.extend_from_slice(qid);
        }
        self.buf.extend_from_slice(line);
        self.buf.push(b'\n');
    }

    /// Detach the closed groups for dispatch, keeping the open group's bytes at
    /// the front of a fresh `buf`. Returns `(bytes, ranges)`; `ranges` index the
    /// returned `bytes`.
    fn take_batch(&mut self) -> (Vec<u8>, Vec<(usize, usize)>) {
        let open = self.buf[self.open_start..].to_vec();
        self.buf.truncate(self.open_start);
        let bytes = std::mem::replace(&mut self.buf, open);
        let ranges = std::mem::take(&mut self.closed);
        self.open_start = 0;
        (bytes, ranges)
    }

    /// At EOF: close the open group and return everything remaining.
    fn finish(&mut self) -> (Vec<u8>, Vec<(usize, usize)>) {
        if self.have_open {
            let end = self.buf.len();
            self.closed.push((self.open_start, end));
            self.have_open = false;
        }
        let bytes = std::mem::take(&mut self.buf);
        let ranges = std::mem::take(&mut self.closed);
        (bytes, ranges)
    }
}

/// Smooth one batch of closed groups in parallel and write the output.
///
/// Worker buffers are collected and emitted in input order. Because batches are
/// processed strictly in sequence, the output is globally order-preserving and
/// byte-identical regardless of `n_threads`.
fn dispatch_batch<W: Write>(
    pool: &rayon::ThreadPool,
    bytes: &[u8],
    ranges: &[(usize, usize)],
    cfg: &SmoothCfg,
    uses_names: bool,
    writer: &mut W,
    stats: &mut SmoothStats,
) {
    use rayon::prelude::*;
    if ranges.is_empty() {
        return;
    }
    let results: Vec<(SmoothStats, Vec<u8>)> = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(s, e)| {
                smooth_group(
                    &bytes[s..e], cfg.tree, cfg.names, cfg.name_to_id, cfg.root_id,
                    uses_names, cfg.max_gap,
                )
            })
            .collect()
    });
    for (st, out) in &results {
        writer.write_all(out).expect("write error");
        stats.merge(st);
    }
}

/// Parallel smoothing pipeline: smooth each query (sequence) on its own thread.
///
/// **Streaming / bounded memory.** The input is consumed in fixed-size blocks
/// and grouped into per-query byte ranges on the fly; complete groups are
/// dispatched to the worker pool in bounded batches, and each batch's output is
/// written and reclaimed before more input is read. Peak memory is therefore
/// `O(batch + largest single query group)` — independent of the total input
/// size — mirroring the streaming property of the Python smoother. It never
/// loads the whole file. A single query group is held whole while it is smoothed
/// (the algorithm needs the full sequence), so an assembly with one enormous
/// chromosome per group peaks at that one group's size, exactly as the Python
/// streamer does.
///
/// Each worker **parses, smooths, and formats its own slice**. The main thread
/// only splits the byte stream into groups (locate newlines and the first tab
/// of each line — no integer parsing or feature resolution), which keeps the two
/// dominant serial costs (parsing and output formatting of tens of millions of
/// lines) inside the parallel region. The earlier version left both on the main
/// thread and parallelised only the smoothing compute, capping speedup near
/// ~1.4x regardless of thread count.
///
/// Output is written in input order and is byte-identical regardless of
/// `n_threads`, so `n_threads == 1` is the canonical reference result. (An
/// order-agnostic mode was measured to give no speedup even on long reads with a
/// 3000x straggler ratio — rayon work-stealing over large batches already
/// balances the load — so it was removed in favour of this single path.)
///
/// `n_threads` controls the worker parallelism (a local rayon pool); it is
/// clamped to at least 1. `n_threads == 1` is the low-memory single-threaded
/// streaming path — there is deliberately no separate single-threaded function
/// to keep in sync, only this one entry point.
pub fn run_smooth(
    input: impl Read,
    output: impl Write,
    tree: &LcaTree,
    names: &[String],
    root_id: usize,
    max_gap: u64,
    n_threads: usize,
) -> SmoothStats {
    let n_threads = n_threads.max(1);
    // Read granularity and batch-dispatch thresholds. A batch closes at whichever
    // limit trips first, keeping peak memory near max_bytes + largest-single-group
    // + output regardless of total input size (the streaming property).
    //
    //  * Reads / many-small-groups: the group cap trips first at a few MiB per
    //    batch, so memory stays tiny while every batch still fills all threads.
    //  * Assemblies / few-huge-groups: the byte cap trips instead; it must be a
    //    healthy multiple of a single group so a batch holds enough chromosomes
    //    to keep the pool busy (a per-group cap would dispatch one chromosome at
    //    a time and serialize them). 512 MiB fits several CHM13 chromosomes yet
    //    stays a constant, input-size-independent bound.
    const READ_BLOCK: usize = 4 * 1024 * 1024; // 4 MiB
    const BATCH_MAX_GROUPS: usize = 65536;
    const BATCH_MAX_BYTES: usize = 512 * 1024 * 1024; // 512 MiB

    let name_to_id: std::collections::HashMap<String, usize> = names
        .iter()
        .enumerate()
        .map(|(id, name)| (name.clone(), id))
        .collect();

    let cfg = SmoothCfg {
        tree,
        names,
        name_to_id: &name_to_id,
        root_id,
        max_gap,
        max_groups: BATCH_MAX_GROUPS,
        max_bytes: BATCH_MAX_BYTES,
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .expect("failed to build rayon thread pool");

    let mut input = input;
    let mut writer = BufWriter::new(output);
    let mut stats = SmoothStats::default();
    let mut grouper = Grouper::default();
    let mut uses_names = false;
    let mut first_content = true;

    // Process one raw line (may span block boundaries): trim, drop blanks,
    // capture/emit the header once, otherwise feed the grouper and dispatch a
    // batch when the thresholds trip.
    let mut handle_line = |raw: &[u8],
                           grouper: &mut Grouper,
                           writer: &mut BufWriter<_>,
                           stats: &mut SmoothStats,
                           uses_names: &mut bool,
                           first_content: &mut bool| {
        let line = raw.trim_ascii();
        if line.is_empty() {
            return;
        }
        if *first_content {
            *first_content = false;
            if line.starts_with(b"query_rank") || line.starts_with(b"query_name") {
                *uses_names = line.windows(b"label_name".len()).any(|w| w == b"label_name");
                writer.write_all(line).expect("write error");
                writer.write_all(b"\n").expect("write error");
                return;
            }
            // Not a header — fall through and treat this line as data.
        }
        grouper.push_line(line);
        if grouper.closed.len() >= cfg.max_groups || grouper.open_start >= cfg.max_bytes {
            let (bytes, ranges) = grouper.take_batch();
            dispatch_batch(&pool, &bytes, &ranges, &cfg, *uses_names, writer, stats);
        }
    };

    let mut line: Vec<u8> = Vec::new();
    let mut block = vec![0u8; READ_BLOCK];
    loop {
        let nread = input.read(&mut block).expect("IO error reading input");
        if nread == 0 {
            break;
        }
        let mut rest = &block[..nread];
        while let Some(pos) = rest.iter().position(|&b| b == b'\n') {
            if line.is_empty() {
                handle_line(&rest[..pos], &mut grouper, &mut writer, &mut stats,
                            &mut uses_names, &mut first_content);
            } else {
                line.extend_from_slice(&rest[..pos]);
                handle_line(&line, &mut grouper, &mut writer, &mut stats,
                            &mut uses_names, &mut first_content);
                line.clear();
            }
            rest = &rest[pos + 1..];
        }
        line.extend_from_slice(rest);
    }
    // Trailing line with no final newline.
    if !line.is_empty() {
        handle_line(&line, &mut grouper, &mut writer, &mut stats,
                    &mut uses_names, &mut first_content);
    }

    // Flush whatever remains.
    let (bytes, ranges) = grouper.finish();
    dispatch_batch(&pool, &bytes, &ranges, &cfg, uses_names, &mut writer, &mut stats);

    writer.flush().expect("flush error");
    stats
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lca_tree::LcaTree;

    /// Build a small tree:   root(3) → A(2) → { B(0), C(1) }
    fn cousin_tree() -> LcaTree {
        // Edges: B→A, C→A, A→root
        LcaTree::new(4, vec![(0, 2), (1, 2), (2, 3)]).unwrap()
    }

    fn iv(feature: usize, start: u64, end: u64) -> Interval {
        Interval { start, end, feature, originally_none: feature == 3 }
    }

    /// Canonical case: B, root, C  →  B, A, C
    /// (a too-general interior interval is promoted to LCA of its specific flankers)
    #[test]
    fn promotes_general_interior_to_lca() {
        let tree = cousin_tree();
        let root = tree.root(); // node 3
        let a = 2usize;
        let b = 0usize;
        let c = 1usize;

        let mut intervals = vec![
            iv(b, 0, 100),
            iv(root, 100, 200),
            iv(c, 200, 300),
        ];

        let reassigned = smooth_intervals(&mut intervals, &tree, 1000);
        assert_eq!(reassigned, 1, "expected exactly one reassignment");
        assert_eq!(intervals[0].feature, b,    "left anchor unchanged");
        assert_eq!(intervals[1].feature, a,    "interior promoted to LCA(B,C) = A");
        assert_eq!(intervals[2].feature, c,    "right anchor unchanged");
    }

    /// Already-at-LCA interior should not be touched.
    #[test]
    fn no_change_when_interior_already_at_lca() {
        let tree = cousin_tree();
        let a = 2usize;
        let b = 0usize;
        let c = 1usize;

        let mut intervals = vec![
            iv(b, 0, 100),
            iv(a, 100, 200),
            iv(c, 200, 300),
        ];

        let reassigned = smooth_intervals(&mut intervals, &tree, 1000);
        assert_eq!(reassigned, 0, "nothing to promote when interior is already LCA");
        assert_eq!(intervals[1].feature, a);
    }

    /// Multi-threaded output must be byte-identical to the single-threaded
    /// (`n_threads == 1`) reference across multiple queries (some needing
    /// promotion, one with a `none` run, one single-interval).
    #[test]
    fn parallel_matches_sequential() {
        let tree = cousin_tree(); // B(0), C(1) → A(2) → root(3)
        let root = tree.root();
        let names = vec![
            "B".to_string(),
            "C".to_string(),
            "A".to_string(),
            "root".to_string(),
        ];
        let input = "query_name\tfrom_kmer\tto_kmer\tlabel_name\n\
                     seq1\t0\t100\tB\n\
                     seq1\t100\t200\troot\n\
                     seq1\t200\t300\tC\n\
                     seq2\t0\t50\tnone\n\
                     seq2\t50\t150\tB\n\
                     seq2\t150\t250\tC\n\
                     seq3\t0\t100\tA\n";

        let mut out_seq: Vec<u8> = Vec::new();
        let s1 = run_smooth(input.as_bytes(), &mut out_seq, &tree, &names, root, 1000, 1);

        for nt in [1usize, 2, 4] {
            let mut out_par: Vec<u8> = Vec::new();
            let s2 =
                run_smooth(input.as_bytes(), &mut out_par, &tree, &names, root, 1000, nt);
            assert_eq!(out_seq, out_par, "parallel output differs at n_threads={nt}");
            assert_eq!(s1.reads_processed, s2.reads_processed);
            assert_eq!(s1.intervals_in, s2.intervals_in);
            assert_eq!(s1.intervals_smoothed, s2.intervals_smoothed);
            assert_eq!(s1.intervals_out, s2.intervals_out);
        }
    }

    /// Force the streaming path across several batch dispatches: with far more
    /// query groups than `BATCH_MAX_GROUPS` (65536), `take_batch` fires mid-stream
    /// multiple times. Output must still be byte-identical to the sequential run,
    /// proving batch boundaries preserve global order.
    #[test]
    fn streaming_batches_preserve_order() {
        let tree = cousin_tree();
        let root = tree.root();
        let names = vec![
            "B".to_string(),
            "C".to_string(),
            "A".to_string(),
            "root".to_string(),
        ];
        // 140000 sequences (> 2 * 65536) each of the promote pattern B, root, C.
        let mut input = String::from("query_name\tfrom_kmer\tto_kmer\tlabel_name\n");
        for i in 0..140000 {
            input.push_str(&format!(
                "seq{i}\t0\t100\tB\nseq{i}\t100\t200\troot\nseq{i}\t200\t300\tC\n"
            ));
        }

        let mut out_seq: Vec<u8> = Vec::new();
        let s1 = run_smooth(input.as_bytes(), &mut out_seq, &tree, &names, root, 1000, 1);

        let mut out_ord: Vec<u8> = Vec::new();
        let s2 =
            run_smooth(input.as_bytes(), &mut out_ord, &tree, &names, root, 1000, 4);
        assert_eq!(out_seq, out_ord, "ordered streaming output differs across batches");
        assert_eq!(s1.reads_processed, s2.reads_processed);
        assert_eq!(s1.reads_processed, 140000);
        assert_eq!(s1.intervals_smoothed, s2.intervals_smoothed);
    }
}
