use std::io::{self, Read, Write};

/// LCA structure with O(n log n) space and O(1) queries.
///
/// Strategy: reduce LCA to Range Minimum Query on the Euler tour.
///   1. DFS the tree, recording each node every time it is visited → Euler tour of length 2n-1.
///   2. `first[v]` = first position of `v` in the tour.
///   3. `LCA(u,v)` = node with minimum depth in `euler[first[u]..=first[v]]`.
///   4. Answer RMQ in O(1) with a sparse table: for each level k store, for every position i,
///      the index of the minimum-depth node in the window `[i, i + 2^k)`.
///      Overlapping windows make queries O(1); building takes O(n log n) time and space.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct LcaSupport {
    n: usize,
    depth: Vec<usize>,
    /// Euler tour: 2n-1 node indices.
    euler: Vec<usize>,
    /// first[v] = first position of node v in `euler`.
    first: Vec<usize>,
    /// Sparse table stored flat: `answer_table[k * m + i]`, where m is the number of nodes
    /// on the Euler tour, is the index into `euler` of the node with minimum depth in the 
    /// window `euler[i .. i + 2^k]` (exclusive endpoint).
    answer_table: Vec<usize>,
}

impl LcaSupport {
    /// Build the structure from `n` nodes and `edges`, where each edge `(child, parent)`
    /// points toward the root. Returns an error if the input is invalid.
    pub fn new(n: usize, edges: Vec<(usize, usize)>) -> Result<Self, String> {
        if n == 0 {
            return Err("Tree must have at least one node".to_string());
        }
        if edges.len() != n - 1 {
            return Err(format!(
                "Expected {} edges for a tree with {} nodes, got {}",
                n - 1,
                n,
                edges.len()
            ));
        }

        let mut parent = vec![usize::MAX; n];
        let mut children: Vec<Vec<usize>> = vec![vec![]; n];

        for &(child, par) in &edges {
            if child >= n || par >= n {
                return Err(format!(
                    "Edge ({child}, {par}) contains out-of-range node index (n={n})"
                ));
            }
            if child == par {
                return Err(format!("Self-loop at node {child}"));
            }
            if parent[child] != usize::MAX {
                return Err(format!("Node {child} has multiple parents"));
            }
            parent[child] = par;
            children[par].push(child);
        }

        let root = (0..n)
            .find(|&i| parent[i] == usize::MAX)
            .ok_or_else(|| "No root found (cycle involving all nodes)".to_string())?;

        // BFS to compute depths and verify connectivity.
        let mut depth = vec![usize::MAX; n];
        depth[root] = 0;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(root);
        let mut visited = 0usize;
        while let Some(node) = queue.pop_front() {
            visited += 1;
            for &child in &children[node] {
                if depth[child] != usize::MAX {
                    return Err(format!("Cycle detected at node {child}"));
                }
                depth[child] = depth[node] + 1;
                queue.push_back(child);
            }
        }
        if visited != n {
            return Err("Tree is disconnected: not all nodes are reachable from root".to_string());
        }

        // Build Euler tour (iterative DFS to avoid stack overflow on deep trees).
        // Each internal visit pushes the node; each return from a child re-pushes the parent.
        // Result: 2n-1 entries.
        let m = 2 * n - 1;
        let mut euler = Vec::with_capacity(m);
        let mut first = vec![0usize; n];
        {
            // Stack entries: (node, next_child_index)
            let mut stack: Vec<(usize, usize)> = Vec::with_capacity(n);
            stack.push((root, 0));
            first[root] = 0;
            euler.push(root);

            while let Some((node, ci)) = stack.last_mut() {
                let node = *node;
                if *ci < children[node].len() {
                    let child = children[node][*ci];
                    *ci += 1;
                    first[child] = euler.len();
                    euler.push(child);
                    stack.push((child, 0));
                } else {
                    stack.pop();
                    // Re-visit parent after returning from this subtree.
                    if let Some(&(parent_node, _)) = stack.last() {
                        euler.push(parent_node);
                    }
                }
            }
        }
        debug_assert_eq!(euler.len(), m);

        // Build answer_table table.
        // answer_table[k * m + i] = index into `euler` of the minimum-depth node in [i, i + 2^k].
        let levels = floor_log2(m) + 1;
        let mut answer_table = vec![0usize; levels * m];

        // Level 0: each window is a single element.
        for i in 0..m {
            answer_table[i] = i;
        }
        // Level k: combine two windows of size 2^(k-1).
        for k in 1..levels {
            let half = 1 << (k - 1);
            for i in 0..m {
                let j = i + half;
                let left = answer_table[(k - 1) * m + i];
                if j >= m {
                    // Right window would be out of bounds; just carry left.
                    answer_table[k * m + i] = left;
                } else {
                    let right = answer_table[(k - 1) * m + j];
                    answer_table[k * m + i] =
                        if depth[euler[left]] <= depth[euler[right]] { left } else { right };
                }
            }
        }

        Ok(LcaSupport { n, depth, euler, first, answer_table })
    }

    /// Returns the depth of `node` in the tree.
    #[inline]
    pub fn depth(&self, node: usize) -> usize {
        self.depth[node]
    }

    /// Returns the lowest common ancestor of nodes `a` and `b`.
    #[inline]
    pub fn lca(&self, a: usize, b: usize) -> usize {
        assert!(a < self.n && b < self.n, "Node index out of bounds");
        if a == b { return a; } // Save some unnecessary work

        // The LCA of a and b is the minimum-depth node in the Euler tour between
        // the first occurrences of a and b.

        // Get the first occurrences of a and b on the Euler tour.
        let (l, r) = {
            let fa = self.first[a];
            let fb = self.first[b];
            if fa <= fb { (fa, fb) } else { (fb, fa) }
        };

        // Look up the precomputed minima for the two overlapping windows of size 2^k
        // that together cover [l, r] exactly.
        let len = r - l + 1;
        let k = floor_log2(len);
        let m = self.euler.len();
        let left  = self.answer_table[k * m + l];
        let right = self.answer_table[k * m + r + 1 - (1 << k)];

        // Return the node at the position with the smaller depth (the "argmin").
        let best = if self.depth[self.euler[left]] <= self.depth[self.euler[right]] {
            left
        } else {
            right
        };

        self.euler[best]
    }

    // --- Serialization ---

    pub fn serialize<W: Write>(&self, w: &mut W) -> io::Result<()> {
        write_u64(w, self.n as u64)?;
        write_usize_slice(w, &self.depth)?;
        write_usize_slice(w, &self.euler)?;
        write_usize_slice(w, &self.first)?;
        write_usize_slice(w, &self.answer_table)?;
        Ok(())
    }

    pub fn load<R: Read>(r: &mut R) -> io::Result<Self> {
        let n = read_u64(r)? as usize;
        let depth = read_usize_vec(r)?;
        let euler = read_usize_vec(r)?;
        let first = read_usize_vec(r)?;
        let answer_table = read_usize_vec(r)?;
        Ok(LcaSupport { n, depth, euler, first, answer_table })
    }
}

#[inline(always)]
fn floor_log2(n: usize) -> usize {
    (usize::BITS - n.leading_zeros() - 1) as usize
}

// --- Binary I/O helpers ---

fn write_u64<W: Write>(w: &mut W, v: u64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_u64<R: Read>(r: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn write_usize_slice<W: Write>(w: &mut W, v: &[usize]) -> io::Result<()> {
    write_u64(w, v.len() as u64)?;
    for &x in v {
        write_u64(w, x as u64)?;
    }
    Ok(())
}

fn read_usize_vec<R: Read>(r: &mut R) -> io::Result<Vec<usize>> {
    let len = read_u64(r)? as usize;
    let mut v = Vec::with_capacity(len);
    for _ in 0..len {
        v.push(read_u64(r)? as usize);
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn binary_tree() -> LcaSupport {
        //         6        <- root (internal)
        //        / \
        //       4   5      <- internal
        //      / \ / \
        //     0  1 2  3   <- leaves
        LcaSupport::new(7, vec![(0,4),(1,4),(2,5),(3,5),(4,6),(5,6)]).unwrap()
    }

    fn path_tree(n: usize) -> LcaSupport {
        // Chain: 0 -> 1 -> 2 -> ... -> n-1 (root)
        let edges = (0..n-1).map(|i| (i, i+1)).collect();
        LcaSupport::new(n, edges).unwrap()
    }

    #[test]
    fn single_node() {
        let t = LcaSupport::new(1, vec![]).unwrap();
        assert_eq!(t.lca(0, 0), 0);
    }

    #[test]
    fn star_lca() {
        // Leaves: 0,1,2,3  Internal/root: 4
        let t = LcaSupport::new(5, vec![(0,4),(1,4),(2,4),(3,4)]).unwrap();
        assert_eq!(t.lca(0, 1), 4);
        assert_eq!(t.lca(2, 3), 4);
        assert_eq!(t.lca(0, 4), 4);
        assert_eq!(t.lca(1, 1), 1);
    }

    #[test]
    fn binary_tree_lca() {
        let t = binary_tree();
        assert_eq!(t.lca(0, 1), 4);
        assert_eq!(t.lca(2, 3), 5);
        assert_eq!(t.lca(0, 2), 6);
        assert_eq!(t.lca(1, 3), 6);
        assert_eq!(t.lca(0, 4), 4);
        assert_eq!(t.lca(4, 5), 6);
        assert_eq!(t.lca(6, 0), 6);
        assert_eq!(t.lca(3, 3), 3);
    }

    #[test]
    fn path_tree_lca() {
        // Chain: 0 -> 1 -> 2 -> 3 -> 4 (root)
        let t = path_tree(5);
        assert_eq!(t.lca(0, 1), 1);
        assert_eq!(t.lca(0, 4), 4);
        assert_eq!(t.lca(1, 3), 3);
        assert_eq!(t.lca(2, 2), 2);
    }

    #[test]
    fn lca_is_symmetric() {
        let t = binary_tree();
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(t.lca(i, j), t.lca(j, i), "LCA({i},{j}) != LCA({j},{i})");
            }
        }
    }

    #[test]
    fn lca_with_self_is_self() {
        let t = binary_tree();
        for i in 0..7 {
            assert_eq!(t.lca(i, i), i);
        }
    }

    #[test]
    fn serialize_roundtrip() {
        let t = binary_tree();
        let mut buf = Vec::new();
        t.serialize(&mut buf).unwrap();
        let t2 = LcaSupport::load(&mut buf.as_slice()).unwrap();

        assert_eq!(t2.n, 7);
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(t2.lca(i, j), t.lca(i, j), "LCA({i},{j}) mismatch after roundtrip");
            }
        }
    }

    #[test]
    fn serialize_roundtrip_path() {
        let t = path_tree(10);
        let mut buf = Vec::new();
        t.serialize(&mut buf).unwrap();
        let t2 = LcaSupport::load(&mut buf.as_slice()).unwrap();

        assert_eq!(t2.n, 10);
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(t2.lca(i, j), t.lca(i, j));
            }
        }
    }

    // --- Validation error cases (mirrors lca_tree.rs) ---

    #[test]
    fn error_wrong_edge_count() {
        let err = LcaSupport::new(4, vec![(1, 3), (2, 3)]).unwrap_err();
        assert!(err.contains("edges"), "{err}");
    }

    #[test]
    fn error_self_loop() {
        let err = LcaSupport::new(3, vec![(0, 2), (2, 2)]).unwrap_err();
        assert!(err.contains("Self-loop"), "{err}");
    }

    #[test]
    fn error_multiple_parents() {
        let err = LcaSupport::new(4, vec![(0, 2), (0, 3), (2, 3)]).unwrap_err();
        assert!(err.contains("multiple parents"), "{err}");
    }

    #[test]
    fn error_disconnected() {
        let err = LcaSupport::new(4, vec![(1, 2), (2, 3), (3, 1)]).unwrap_err();
        assert!(err.contains("disconnected"), "{err}");
    }

    #[test]
    fn error_out_of_range_node() {
        let err = LcaSupport::new(3, vec![(0, 2), (5, 2)]).unwrap_err();
        assert!(err.contains("out-of-range"), "{err}");
    }

    #[test]
    fn error_zero_nodes() {
        let err = LcaSupport::new(0, vec![]).unwrap_err();
        assert!(err.contains("at least one"), "{err}");
    }

    /// Cross-check: LcaSupport must agree with the naive O(n²) LcaTree on a larger random-ish tree.
    #[test]
    fn agrees_with_naive_on_larger_tree() {
        use crate::lca_tree::LcaTree;

        // Build a balanced-ish tree with 15 nodes (4-level complete binary tree).
        // Leaves: 0..7, internals: 8..14, root: 14.
        //           14
        //         /    \
        //       12      13
        //      /  \    /  \
        //     8    9  10   11
        //    / \ / \ / \  / \
        //    0 1 2 3 4 5  6  7
        let edges = vec![
            (0,8),(1,8),(2,9),(3,9),(4,10),(5,10),(6,11),(7,11),
            (8,12),(9,12),(10,13),(11,13),(12,14),(13,14),
        ];
        let fast = LcaSupport::new(15, edges.clone()).unwrap();
        let naive = LcaTree::new(15, edges).unwrap();

        for i in 0..15 {
            for j in 0..15 {
                assert_eq!(
                    fast.lca(i, j), naive.lca(i, j),
                    "Disagreement at LCA({i},{j})"
                );
            }
        }
    }
}
