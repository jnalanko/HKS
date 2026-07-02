//! Hierarchy-aware smoothing for HKS lookup output.
//!
//! Implements Algorithm S1 from the HKS paper: a two-phase scan that identifies
//! windows exhibiting a specific → general → specific pattern in the category
//! hierarchy and reassigns interior intervals to the LCA of the flanking anchors.

use std::collections::HashSet;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

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

/// Flush a completed query: smooth → merge → write.
fn flush_query(
    query_id: &str,
    buf: &mut Vec<Interval>,
    tree: &LcaTree,
    names: &[String],
    root_id: usize,
    uses_names: bool,
    max_gap: u64,
    writer: &mut impl Write,
    stats: &mut SmoothStats,
) {
    let n_in = buf.len() as u64;
    let reassigned = smooth_intervals(buf, tree, max_gap);
    let (merged, eliminated) = merge_intervals(std::mem::take(buf));
    let n_out = merged.len() as u64;

    for iv in &merged {
        let feat_str = format_feature(iv.feature, iv.originally_none, names, root_id, uses_names);
        writeln!(writer, "{}\t{}\t{}\t{}", query_id, iv.start, iv.end, feat_str)
            .expect("write error");
    }

    stats.reads_processed += 1;
    stats.intervals_in += n_in;
    stats.intervals_smoothed += reassigned;
    stats.intervals_merged += eliminated;
    stats.intervals_out += n_out;
}

/// Run the smoothing pipeline on TSV input.
pub fn run_smooth(
    input: impl Read,
    output: impl Write,
    tree: &LcaTree,
    names: &[String],
    root_id: usize,
    max_gap: u64,
) -> SmoothStats {
    // Build name → id lookup
    let name_to_id: std::collections::HashMap<String, usize> = names.iter()
        .enumerate()
        .map(|(id, name)| (name.clone(), id))
        .collect();

    let reader = BufReader::new(input);
    let mut writer = BufWriter::new(output);
    let mut stats = SmoothStats::default();
    let mut buf: Vec<Interval> = Vec::new();
    let mut current_query = String::new();
    let mut uses_names = false; // determined from header
    let mut header_seen = false;

    for line in reader.lines() {
        let line = line.expect("IO error reading input");
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        // Parse header
        if !header_seen {
            if trimmed.starts_with("query_rank") || trimmed.starts_with("query_name") {
                header_seen = true;
                uses_names = trimmed.contains("label_name");
                // Pass header through
                writeln!(writer, "{}", trimmed).expect("write error");
                continue;
            }
            // No header line — treat as data
            header_seen = true;
        }

        let mut cols = trimmed.splitn(4, '\t');
        let query_id = cols.next().expect("missing query column");
        let start: u64 = cols.next()
            .and_then(|s| s.parse().ok())
            .expect("bad start coordinate");
        let end: u64 = cols.next()
            .and_then(|s| s.parse().ok())
            .expect("bad end coordinate");
        let feat_token = cols.next().expect("missing feature column");

        let (feature, originally_none) = resolve_feature(feat_token, &name_to_id, root_id, uses_names);

        // Flush on query boundary
        if query_id != current_query {
            if !buf.is_empty() {
                log::info!("Smoothing {}", current_query);
                flush_query(&current_query, &mut buf, tree, names, root_id, uses_names, max_gap, &mut writer, &mut stats);
            }
            current_query = query_id.to_string();
        }

        buf.push(Interval { start, end, feature, originally_none });
    }

    // Flush final query
    if !buf.is_empty() {
        log::info!("Smoothing {}", current_query);
        flush_query(&current_query, &mut buf, tree, names, root_id, uses_names, max_gap, &mut writer, &mut stats);
    }

    writer.flush().expect("flush error");
    stats
}

/// Parallel smoothing pipeline: smooth each query (sequence) on its own thread.
///
/// Each query's intervals are independent, so smoothing + merging can run in
/// parallel across queries. This function:
///   1. parses the whole input into consecutive per-query interval groups
///      (preserving input order, exactly like [`run_smooth`]'s flush-on-change),
///   2. smooths + merges each group in parallel via rayon, then
///   3. writes the groups back out in the original input order.
///
/// The output is byte-for-byte identical to [`run_smooth`]; only the work is
/// distributed. `n_threads` controls the smoothing parallelism (a local rayon
/// pool is built for it). Memory: unlike the streaming [`run_smooth`], this
/// holds the whole input's intervals in memory at once (still far smaller than
/// the HKS index itself).
pub fn run_smooth_parallel(
    input: impl Read,
    output: impl Write,
    tree: &LcaTree,
    names: &[String],
    root_id: usize,
    max_gap: u64,
    n_threads: usize,
) -> SmoothStats {
    use rayon::prelude::*;

    // Build name → id lookup
    let name_to_id: std::collections::HashMap<String, usize> = names.iter()
        .enumerate()
        .map(|(id, name)| (name.clone(), id))
        .collect();

    let reader = BufReader::new(input);
    let mut writer = BufWriter::new(output);

    // ---- Phase 1: parse into consecutive per-query groups, in order ----
    let mut groups: Vec<(String, Vec<Interval>)> = Vec::new();
    let mut uses_names = false;
    let mut header_seen = false;
    let mut header_line: Option<String> = None;

    for line in reader.lines() {
        let line = line.expect("IO error reading input");
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if !header_seen {
            if trimmed.starts_with("query_rank") || trimmed.starts_with("query_name") {
                header_seen = true;
                uses_names = trimmed.contains("label_name");
                header_line = Some(trimmed.to_string());
                continue;
            }
            header_seen = true;
        }

        let mut cols = trimmed.splitn(4, '\t');
        let query_id = cols.next().expect("missing query column");
        let start: u64 = cols.next()
            .and_then(|s| s.parse().ok())
            .expect("bad start coordinate");
        let end: u64 = cols.next()
            .and_then(|s| s.parse().ok())
            .expect("bad end coordinate");
        let feat_token = cols.next().expect("missing feature column");

        let (feature, originally_none) = resolve_feature(feat_token, &name_to_id, root_id, uses_names);

        // Start a new group when the query id changes (flush-on-change, so the
        // grouping matches run_smooth's byte output for consecutive queries).
        if groups.last().map(|(q, _)| q.as_str()) != Some(query_id) {
            groups.push((query_id.to_string(), Vec::new()));
        }
        groups.last_mut().unwrap().1.push(Interval { start, end, feature, originally_none });
    }

    // Header goes out first, matching run_smooth.
    if let Some(h) = &header_line {
        writeln!(writer, "{}", h).expect("write error");
    }

    // ---- Phase 2: smooth + merge each group in parallel ----
    // Per-group result: (n_in, reassigned, eliminated, merged intervals).
    // Only this parallel map runs in the pool; the (non-Send) input/output
    // handles are touched solely on this (main) thread, before and after.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .expect("failed to build rayon thread pool");
    let results: Vec<(u64, u64, u64, Vec<Interval>)> = pool.install(|| {
        groups
            .par_iter_mut()
            .map(|(_qid, buf)| {
                let n_in = buf.len() as u64;
                let reassigned = smooth_intervals(buf, tree, max_gap);
                let (merged, eliminated) = merge_intervals(std::mem::take(buf));
                (n_in, reassigned, eliminated, merged)
            })
            .collect()
    });

    // ---- Phase 3: write in input order + accumulate stats ----
    let mut stats = SmoothStats::default();
    for ((qid, _), (n_in, reassigned, eliminated, merged)) in groups.iter().zip(results.iter()) {
        for iv in merged {
            let feat_str = format_feature(iv.feature, iv.originally_none, names, root_id, uses_names);
            writeln!(writer, "{}\t{}\t{}\t{}", qid, iv.start, iv.end, feat_str)
                .expect("write error");
        }
        stats.reads_processed += 1;
        stats.intervals_in += n_in;
        stats.intervals_smoothed += reassigned;
        stats.intervals_merged += eliminated;
        stats.intervals_out += merged.len() as u64;
    }

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

    /// The parallel pipeline must produce byte-identical output to the
    /// sequential one across multiple queries (some needing promotion, one
    /// with a `none` run, one single-interval).
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
        let s1 = run_smooth(input.as_bytes(), &mut out_seq, &tree, &names, root, 1000);

        for nt in [1usize, 2, 4] {
            let mut out_par: Vec<u8> = Vec::new();
            let s2 = run_smooth_parallel(input.as_bytes(), &mut out_par, &tree, &names, root, 1000, nt);
            assert_eq!(out_seq, out_par, "parallel output differs at n_threads={nt}");
            assert_eq!(s1.reads_processed, s2.reads_processed);
            assert_eq!(s1.intervals_in, s2.intervals_in);
            assert_eq!(s1.intervals_smoothed, s2.intervals_smoothed);
            assert_eq!(s1.intervals_out, s2.intervals_out);
        }
    }
}
