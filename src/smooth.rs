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

/// Reusable working buffers for [`smooth_intervals`].
///
/// The window scan allocates per *window*, and a whole-chromosome query group
/// contains millions of windows, so allocating these fresh each time dominates
/// allocator traffic in a parallel run. Hoisting them here lets one buffer set
/// serve every window of every group a worker handles: they are cleared, not
/// reallocated, once they reach steady-state capacity.
#[derive(Default)]
pub struct SmoothScratch {
    was_related: Vec<bool>,
    related: Vec<usize>,
    disallowed: HashSet<usize>,
}

/// Smooth a single query's intervals in-place using the hierarchy.
/// Returns the number of feature reassignments made.
///
/// Allocates its own scratch. Use [`smooth_intervals_with`] when smoothing many
/// queries so the buffers are reused across them.
pub fn smooth_intervals(intervals: &mut Vec<Interval>, tree: &LcaTree, max_gap: u64) -> u64 {
    smooth_intervals_with(intervals, tree, max_gap, &mut SmoothScratch::default())
}

/// [`smooth_intervals`], reusing caller-owned buffers.
pub fn smooth_intervals_with(
    intervals: &mut Vec<Interval>,
    tree: &LcaTree,
    max_gap: u64,
    scratch: &mut SmoothScratch,
) -> u64 {
    if intervals.len() < 2 {
        return 0;
    }
    // Split the borrow up front so all three buffers can be held at once.
    let SmoothScratch { was_related, related, disallowed } = scratch;
    let mut total_reassignments = 0u64;
    let n = intervals.len();

    loop {
        let mut changed = false;
        was_related.clear();
        was_related.resize(n, false);
        let mut i = 0;

        while i < n {
            let window_start = i;
            related.clear();
            related.push(i);
            let mut first_unrelated: Option<usize> = None;

            // ------------------------------------------------------------------
            // Ascending phase: scan rightward accepting features that are
            // ancestors of last_rel_feat (i.e. more general / closer to root).
            // Unrelated features on other branches are skipped but their ancestor
            // paths are added to `disallowed` so we stop if we'd cross into them.
            // ------------------------------------------------------------------
            disallowed.clear();
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
            for &idx in related.iter() {
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
pub fn merge_intervals(intervals: &mut Vec<Interval>) -> u64 {
    let n_in = intervals.len();
    // `dedup_by` passes (next, current) and drops `next` when the closure is
    // true, so extending `current` there is exactly the merge. Done in place:
    // this runs on every worker for every group, and the old version allocated a
    // second Vec the size of the input each time.
    intervals.dedup_by(|next, current| {
        if next.feature == current.feature
            && next.start == current.end
            && next.originally_none == current.originally_none
        {
            current.end = next.end;
            true
        } else {
            false
        }
    });
    (n_in - intervals.len()) as u64
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

/// The label that marks a miss in the input, which is also the one we write back
/// out. This is the value of --miss-label if it was given, and otherwise the
/// default of the lookup command, which depends on whether the input has label
/// names or label ids.
fn resolve_miss_label(miss_label: Option<&str>, uses_names: bool) -> &str {
    miss_label.unwrap_or(if uses_names { crate::DEFAULT_MISS_LABEL } else { crate::DEFAULT_MISS_LABEL_ID_MODE })
}

/// Resolve a feature token from the TSV input into a node ID.
/// Returns (node_id, originally_none).
fn resolve_feature(
    token: &str,
    name_to_id: &std::collections::HashMap<String, usize>,
    root_id: usize,
    uses_names: bool,
    miss_label: &str,
) -> (usize, bool) {
    if token == miss_label {
        (root_id, true)
    } else if uses_names {
        let id = name_to_id.get(token)
            .unwrap_or_else(|| panic!("Unknown feature name in input: {:?}", token));
        (*id, false)
    } else {
        // Numeric ID mode
        let id: usize = token.parse()
            .unwrap_or_else(|_| panic!("Cannot parse feature ID: {:?}", token));
        (id, false)
    }
}

/// Append a feature token to the output buffer.
///
/// Writes bytes straight into `out` rather than returning an owned `String`.
/// This is called once per *output interval*, so returning a `String` meant one
/// short-lived heap allocation per output line — tens of millions per run, from
/// every worker at once.
fn write_feature(
    out: &mut Vec<u8>,
    feature: usize,
    originally_none: bool,
    names: &[String],
    root_id: usize,
    uses_names: bool,
    miss_label: &str,
) {
    if originally_none && feature == root_id {
        // Was a miss and smoothing didn't resolve it → keep it as a miss
        out.extend_from_slice(miss_label.as_bytes());
    } else if uses_names {
        out.extend_from_slice(names[feature].as_bytes());
    } else {
        write!(out, "{feature}").expect("formatting error");
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
    cfg: &SmoothCfg,
    uses_names: bool,
    intervals: &mut Vec<Interval>,
    smooth_scratch: &mut SmoothScratch,
    out: &mut Vec<u8>,
) -> SmoothStats {
    let mut query_id: &str = "";
    let mut have_qid = false;
    intervals.clear();
    out.clear();

    let miss_label = resolve_miss_label(cfg.miss_label, uses_names);

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

        let (feature, originally_none) =
            resolve_feature(feat_token, cfg.name_to_id, cfg.root_id, uses_names, miss_label);
        intervals.push(Interval { start, end, feature, originally_none });
    }

    let n_in = intervals.len() as u64;
    let reassigned = smooth_intervals_with(intervals, cfg.tree, cfg.max_gap, smooth_scratch);
    let eliminated = merge_intervals(intervals);
    let n_out = intervals.len() as u64;

    // `out` is a recycled buffer that has already grown to whatever this
    // workload needs, so it is neither allocated nor sized here. Sizing per
    // group is a trap in both directions: reserving the input size mmaps
    // hundreds of MB per worker for a group that collapses to a few intervals,
    // while a fixed reserve is huge overhead for read data, where a group's
    // output is a few dozen bytes and there are millions of them.
    for iv in intervals.iter() {
        write!(out, "{}\t{}\t{}\t", query_id, iv.start, iv.end).expect("formatting error");
        write_feature(out, iv.feature, iv.originally_none, cfg.names, cfg.root_id, uses_names, miss_label);
        out.push(b'\n');
    }

    SmoothStats {
        reads_processed: 1,
        intervals_in: n_in,
        intervals_smoothed: reassigned,
        intervals_merged: eliminated,
        intervals_out: n_out,
    }
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
    miss_label: Option<&'a str>, // None = the lookup command default, see resolve_miss_label
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

thread_local! {
    /// Per-worker parse and smoothing buffers, one set per pool thread for the
    /// whole run.
    ///
    /// Deliberately not `map_init`: rayon builds that closure's value once per
    /// *leaf of the split tree*, not once per thread, and the leaf count climbs
    /// steeply with thread count -- measured at 4 leaves for 1 thread but 3190
    /// for 16 on one short-read input. Each fresh set then has to re-grow its
    /// buffers, which for query groups small enough that a leaf covers only a
    /// couple of dozen of them costs more than the reuse saves. Keying off the
    /// thread bounds the number of buffer sets by the pool size.
    static WORKER_BUFS: std::cell::RefCell<(Vec<Interval>, SmoothScratch)> =
        std::cell::RefCell::new((Vec::new(), SmoothScratch::default()));
}

/// Smooth one batch of closed groups in parallel and write the output.
///
/// `outbufs` is the run's recycled output-buffer pool, owned by `run_smooth` so
/// it survives across batches; it grows only to the largest group count any one
/// batch needs, and each buffer keeps the capacity it reached.
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
    outbufs: &mut Vec<Vec<u8>>,
) {
    use rayon::prelude::*;
    if ranges.is_empty() {
        return;
    }
    if outbufs.len() < ranges.len() {
        outbufs.resize_with(ranges.len(), Vec::new);
    }
    // One output buffer per group, taken from the run's pool rather than
    // allocated. `Zip` is an indexed parallel iterator, so results still come
    // back in input order and the output stays byte-identical at any thread
    // count.
    let results: Vec<SmoothStats> = pool.install(|| {
        ranges
            .par_iter()
            .zip(outbufs[..ranges.len()].par_iter_mut())
            .map(|(&(s, e), out)| {
                WORKER_BUFS.with(|cell| {
                    let (intervals, smooth_scratch) = &mut *cell.borrow_mut();
                    smooth_group(&bytes[s..e], cfg, uses_names, intervals, smooth_scratch, out)
                })
            })
            .collect()
    });
    for (st, out) in results.iter().zip(outbufs.iter()) {
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
///
/// With `no_header` the first line is data rather than a header (e.g. a BED
/// file), and nothing is written to the output before the first record. The
/// label column format then comes from `report_label_ids` instead of from the
/// header.
pub fn run_smooth(
    input: impl Read,
    output: impl Write,
    tree: &LcaTree,
    names: &[String],
    root_id: usize,
    max_gap: u64,
    miss_label: Option<&str>,
    no_header: bool,
    report_label_ids: bool,
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
        miss_label,
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
    // Recycled across every batch; see dispatch_batch.
    let mut outbufs: Vec<Vec<u8>> = Vec::new();
    // A header line, if there is one, states whether the label column has names
    // or ids. With --no-header there is no line to sniff, so the caller says
    // which it is. Without a header line and without --no-header we keep the old
    // fallback of treating the input as ids.
    let mut uses_names = no_header && !report_label_ids;
    let mut first_content = !no_header;

    // Process one raw line (may span block boundaries): trim, drop blanks,
    // capture/emit the header once, otherwise feed the grouper and dispatch a
    // batch when the thresholds trip.
    let mut handle_line = |raw: &[u8],
                           grouper: &mut Grouper,
                           writer: &mut BufWriter<_>,
                           stats: &mut SmoothStats,
                           uses_names: &mut bool,
                           first_content: &mut bool,
                           outbufs: &mut Vec<Vec<u8>>| {
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
            dispatch_batch(&pool, &bytes, &ranges, &cfg, *uses_names, writer, stats, outbufs);
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
                            &mut uses_names, &mut first_content, &mut outbufs);
            } else {
                line.extend_from_slice(&rest[..pos]);
                handle_line(&line, &mut grouper, &mut writer, &mut stats,
                            &mut uses_names, &mut first_content, &mut outbufs);
                line.clear();
            }
            rest = &rest[pos + 1..];
        }
        line.extend_from_slice(rest);
    }
    // Trailing line with no final newline.
    if !line.is_empty() {
        handle_line(&line, &mut grouper, &mut writer, &mut stats,
                    &mut uses_names, &mut first_content, &mut outbufs);
    }

    // Flush whatever remains.
    let (bytes, ranges) = grouper.finish();
    dispatch_batch(&pool, &bytes, &ranges, &cfg, uses_names, &mut writer, &mut stats, &mut outbufs);

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
        let s1 = run_smooth(input.as_bytes(), &mut out_seq, &tree, &names, root, 1000, None, false, false, 1);

        for nt in [1usize, 2, 4] {
            let mut out_par: Vec<u8> = Vec::new();
            let s2 =
                run_smooth(input.as_bytes(), &mut out_par, &tree, &names, root, 1000, None, false, false, nt);
            assert_eq!(out_seq, out_par, "parallel output differs at n_threads={nt}");
            assert_eq!(s1.reads_processed, s2.reads_processed);
            assert_eq!(s1.intervals_in, s2.intervals_in);
            assert_eq!(s1.intervals_smoothed, s2.intervals_smoothed);
            assert_eq!(s1.intervals_out, s2.intervals_out);
        }
    }

    /// Headerless input (e.g. a BED file) must smooth to the same records as the
    /// same input with a header, minus the header line itself.
    #[test]
    fn headerless_input_matches_headered() {
        let tree = cousin_tree();
        let root = tree.root();
        let names = vec![
            "B".to_string(),
            "C".to_string(),
            "A".to_string(),
            "root".to_string(),
        ];
        let records = "seq1\t0\t100\tB\n\
                       seq1\t100\t200\troot\n\
                       seq1\t200\t300\tC\n\
                       seq2\t0\t50\tnone\n\
                       seq2\t50\t150\tB\n";
        let with_header = format!("query_name\tfrom_kmer\tto_kmer\tlabel_name\n{records}");

        let mut out_headered: Vec<u8> = Vec::new();
        run_smooth(with_header.as_bytes(), &mut out_headered, &tree, &names, root, 1000, None, false, false, 1);

        for nt in [1usize, 2] {
            let mut out_headerless: Vec<u8> = Vec::new();
            run_smooth(records.as_bytes(), &mut out_headerless, &tree, &names, root, 1000, None, true, false, nt);
            let headered_records = &out_headered[out_headered.iter().position(|&b| b == b'\n').unwrap() + 1..];
            assert_eq!(headered_records, out_headerless, "headerless output differs at n_threads={nt}");
        }
    }

    /// Headerless input in label id mode: the label column is parsed as node ids,
    /// and an unresolved miss keeps the id-mode miss label.
    #[test]
    fn headerless_input_with_label_ids() {
        let tree = cousin_tree(); // B(0), C(1) → A(2) → root(3)
        let root = tree.root();
        let names = vec![
            "B".to_string(),
            "C".to_string(),
            "A".to_string(),
            "root".to_string(),
        ];
        let records = "seq1\t0\t100\t0\n\
                       seq1\t100\t200\t3\n\
                       seq1\t200\t300\t1\n\
                       seq2\t0\t50\t-\n";

        let mut out: Vec<u8> = Vec::new();
        run_smooth(records.as_bytes(), &mut out, &tree, &names, root, 1000, None, true, true, 1);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "seq1\t0\t100\t0\n\
             seq1\t100\t200\t2\n\
             seq1\t200\t300\t1\n\
             seq2\t0\t50\t-\n",
            "interior root(3) promoted to LCA(B,C) = A(2), miss kept as '-'"
        );
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
        let s1 = run_smooth(input.as_bytes(), &mut out_seq, &tree, &names, root, 1000, None, false, false, 1);

        let mut out_ord: Vec<u8> = Vec::new();
        let s2 =
            run_smooth(input.as_bytes(), &mut out_ord, &tree, &names, root, 1000, None, false, false, 4);
        assert_eq!(out_seq, out_ord, "ordered streaming output differs across batches");
        assert_eq!(s1.reads_processed, s2.reads_processed);
        assert_eq!(s1.reads_processed, 140000);
        assert_eq!(s1.intervals_smoothed, s2.intervals_smoothed);
    }
}
