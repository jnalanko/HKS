//! Hierarchy-aware smoothing for HKS lookup output.
//!
//! Implements Algorithm S1 from the HKS paper: a two-phase scan that identifies
//! windows exhibiting a specific → general → specific pattern in the category
//! hierarchy and reassigns interior intervals to the LCA of the flanking anchors.

use std::collections::HashSet;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use crate::lca_tree::LcaTree;
use crate::single_colored_kmers::ColorHierarchy;

// ---------------------------------------------------------------------------
// Interval representation
// ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct Interval {
    pub query_id: String,
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
            // Ancestor phase: scan rightward accepting intervals whose features
            // are ancestors of the left anchor's feature (more general).
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

                if tree.is_ancestor(last_rel_feat, nf) {
                    // nf IS an ancestor of last_rel_feat → more general → extend
                    if disallowed.contains(&nf) {
                        break;
                    }
                    last_rel_feat = nf;
                    last_rel_end = intervals[j + 1].end;
                    last_rel_idx = j + 1;
                    related.push(j + 1);
                    j += 1;
                } else if !tree.is_ancestor(nf, last_rel_feat) {
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
            // Descendant phase: continue rightward from the turning point,
            // accepting intervals whose features are descendants of the turning
            // point's feature (more specific).
            // ------------------------------------------------------------------
            let peak_feat = intervals[window_end].feature;
            let mut k = window_end;

            while k + 1 < n {
                let nf = intervals[k + 1].feature;
                let ns = intervals[k + 1].start;

                if ns > last_rel_end + max_gap {
                    break;
                }

                if tree.is_ancestor(nf, last_rel_feat) {
                    // last_rel_feat IS ancestor of nf → nf is more specific → extend
                    last_rel_feat = nf;
                    last_rel_end = intervals[k + 1].end;
                    last_rel_idx = k + 1;
                    related.push(k + 1);
                    k += 1;
                } else if !tree.is_ancestor(last_rel_feat, nf) {
                    // Neither is ancestor → unrelated
                    if tree.is_ancestor(nf, peak_feat) {
                        // nf is a descendant of peak → would restart a new window → stop
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
                    // is_ancestor(lca, orig) means orig IS an ancestor of lca,
                    // i.e. orig is more general than lca → replace with lca
                    if tree.is_ancestor(lca, orig) && orig != lca {
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
    buf: &mut Vec<Interval>,
    hierarchy: &ColorHierarchy,
    uses_names: bool,
    max_gap: u64,
    writer: &mut impl Write,
    stats: &mut SmoothStats,
) {
    let n_in = buf.len() as u64;
    let reassigned = smooth_intervals(buf, hierarchy.tree(), max_gap);
    let (merged, eliminated) = merge_intervals(std::mem::take(buf));
    let n_out = merged.len() as u64;

    for iv in &merged {
        let feat_str = format_feature(iv.feature, iv.originally_none, hierarchy.names(), hierarchy.root(), uses_names);
        writeln!(writer, "{}\t{}\t{}\t{}", iv.query_id, iv.start, iv.end, feat_str)
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
    hierarchy: &ColorHierarchy,
    max_gap: u64,
) -> SmoothStats {
    let names = hierarchy.names();
    let root_id = hierarchy.root();
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
                uses_names = trimmed.contains("color_name");
                // Pass header through
                writeln!(writer, "{}", trimmed).expect("write error");
                continue;
            }
            // No header line — treat as data
            header_seen = true;
        }

        let mut cols = trimmed.splitn(4, '\t');
        let query_id = cols.next().expect("missing query column").to_string();
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
                flush_query(&mut buf, hierarchy, uses_names, max_gap, &mut writer, &mut stats);
            }
            current_query = query_id.clone();
        }

        buf.push(Interval {
            query_id,
            start,
            end,
            feature,
            originally_none,
        });
    }

    // Flush final query
    if !buf.is_empty() {
        log::info!("Smoothing {}", current_query);
        flush_query(&mut buf, hierarchy, uses_names, max_gap, &mut writer, &mut stats);
    }

    writer.flush().expect("flush error");
    stats
}
