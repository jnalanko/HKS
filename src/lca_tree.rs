use std::io::{self, Read, Write};
use crate::lca_support::LcaSupport;

/// A rooted tree where leaves occupy IDs 0..num_leaves and internal nodes occupy
/// IDs num_leaves..n. The root is an internal node (or node 0 for a single-node tree).
/// `parent[root] == root` (root is its own parent).
///
/// Optionally supports **priority-aware LCA** (`plca`), which gives users
/// more control over how multi-label k-mers are resolved.
///
/// With standard LCA using the chromosome hierarchy, a k-mer found in both
/// "chr1" and "chrM" would be labeled as their common ancestor (e.g. "root"),
/// losing specificity. Using priorities, the user can declare that
/// "mitochondrial" is more interesting than "autosomal" for their use case.
/// A k-mer shared between chr1 and chrM would then be assigned to chrM
/// instead of root.
///
/// For `plca`, each node in the hierarchy is assigned an integer priority
/// (lower value = higher priority). Priorities determine the relative
/// importance of sibling nodes. For example, in the chromosome hierarchy
/// autosomal, mitochondrial, and sex are siblings under the root node. A
/// user could designate their respective priorities as 1, 2, and 3 which
/// would indicate that autosomal has the highest priority and sex has the
/// lowest priority.
///
/// When computing `plca(A, B)`:
/// - If A and B are the same node or one is an ancestor of the other,
///   `plca` returns the same as the standard LCA.
/// - Otherwise, let C_A and C_B be the children of LCA(A, B) that are
///   ancestors of A and B respectively. If C_A and C_B have the same
///   priority (or no priorities are assigned), `plca` returns the same
///   as the standard LCA. If C_A has higher priority (lower number)
///   than C_B, `plca` returns A. If C_B has higher priority, `plca`
///   returns B.
///
/// For any group of siblings, the priorities must be either all the same
/// or all distinct. These conditions ensure that `plca` is associative
/// and safe for use in incremental folds. A mix of same and distinct
/// priorities among siblings breaks associativity; [`set_priorities`]
/// will reject such configurations and suggest either assigning distinct
/// priorities or adding a named grouping node as the parent of the tied
/// siblings.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct LcaTree {
    n: usize,
    root: usize,
    /// parent[i] is the parent of node i; parent[root] == root.
    parent: Vec<usize>,
    lca_support: LcaSupport,
    /// priority[i] is the priority of node i (lower = higher priority).
    /// When empty, `plca` falls back to standard `lca`.
    priority: Vec<usize>,
    /// ancestors_by_depth[i] is a flat list of node i's ancestors from itself up to
    /// the root: `[i, parent(i), grandparent(i), ..., root]`.
    /// Used to answer "child of ancestor A on path to descendant D" in O(1):
    ///   `ancestors_by_depth[D][depth(D) - depth(A) - 1]`
    ancestors_by_depth: Vec<Vec<usize>>,
}

impl LcaTree {
    /// Construct a tree from `n` nodes and `edges`, where each edge `(child, parent)`
    /// points toward the root. Leaves must have IDs 0..num_leaves and internal nodes
    /// must have IDs num_leaves..n. Returns an error if these conditions are not met or
    /// if the edges do not form a valid tree.
    ///
    /// The tree is created without priorities; `plca` will behave identically to `lca`.
    /// Call [`set_priorities`] afterwards to enable priority-aware LCA.
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

        // usize::MAX = "not yet assigned"
        let mut parent = vec![usize::MAX; n];
        let mut children = vec![vec![]; n];

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

        // The root is the unique node with no parent edge.
        let root = (0..n)
            .find(|&i| parent[i] == usize::MAX)
            .ok_or_else(|| "No root found (cycle involving all nodes)".to_string())?;
        parent[root] = root;

        // Validate leaf/internal node ordering: leaves (empty children) must have IDs
        // 0..num_leaves, internal nodes (non-empty children) must have IDs num_leaves..n.
        let num_leaves = children.iter().take_while(|c| c.is_empty()).count();
        for i in num_leaves..n {
            if children[i].is_empty() {
                return Err(format!(
                    "Node {i} is a leaf but has ID >= {num_leaves}: \
                     leaves must occupy IDs 0..{num_leaves}"
                ));
            }
        }

        // BFS from root to compute depths and verify connectivity.
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

        // Build ancestor lists: ancestors_by_depth[v] = [v, parent(v), ..., root].
        let mut ancestors_by_depth = vec![vec![]; n];
        for v in 0..n {
            let mut anc = Vec::with_capacity(depth[v] + 1);
            let mut cur = v;
            loop {
                anc.push(cur);
                if cur == root {
                    break;
                }
                cur = parent[cur];
            }
            ancestors_by_depth[v] = anc;
        }

        let lca_support = LcaSupport::new(n, edges).unwrap();

        Ok(LcaTree { n, root, parent, lca_support, priority: vec![], ancestors_by_depth })
    }

    /// Set per-node priorities for priority-aware LCA.
    ///
    /// `priorities` must have exactly `n` elements, one per node. Lower values
    /// mean higher priority.
    ///
    /// For each group of siblings, the priorities must be either **all the same**
    /// (equivalent to no priority preference — `plca` returns the standard LCA)
    /// or **all distinct** (full priority-aware resolution). A mix of same and
    /// distinct priorities among siblings breaks associativity and is rejected.
    ///
    /// Returns an error if the length is wrong or if any sibling group has a
    /// mixed configuration. The error message lists all violations.
    pub fn set_priorities(&mut self, priorities: Vec<usize>) -> Result<(), String> {
        if priorities.len() != self.n {
            return Err(format!(
                "Expected {} priorities, got {}",
                self.n,
                priorities.len()
            ));
        }

        // Build children lists from parent array for validation.
        let mut children = vec![vec![]; self.n];
        for v in 0..self.n {
            if v != self.root {
                children[self.parent[v]].push(v);
            }
        }

        // For each sibling group, check that priorities are either all the same
        // or all distinct. Collect all violations before reporting.
        let mut violations = Vec::new();
        for v in 0..self.n {
            if children[v].len() < 2 {
                continue;
            }
            let child_priorities: Vec<usize> = children[v].iter().map(|&c| priorities[c]).collect();
            let all_same = child_priorities.iter().all(|&p| p == child_priorities[0]);
            if all_same {
                continue; // All same priority is fine — standard LCA behavior.
            }
            // Not all the same — check if all distinct.
            let mut seen = std::collections::HashSet::new();
            let mut duplicates: std::collections::HashMap<usize, Vec<usize>> =
                std::collections::HashMap::new();
            for &child in &children[v] {
                if !seen.insert(priorities[child]) {
                    duplicates.entry(priorities[child]).or_default().push(child);
                } else {
                    // First occurrence — add it too so the full group is listed.
                    duplicates.entry(priorities[child]).or_default().push(child);
                }
            }
            // Keep only the priorities that have duplicates.
            for (pri, nodes) in &duplicates {
                if nodes.len() > 1 {
                    let node_list: Vec<String> = nodes.iter().map(|n| n.to_string()).collect();
                    violations.push(format!(
                        "  Children of node {v} with duplicate priority {pri}: [{}]. \
                         Either assign them distinct priorities or add a named grouping \
                         node as their parent at priority {pri}.",
                        node_list.join(", ")
                    ));
                }
            }
        }

        if !violations.is_empty() {
            return Err(format!(
                "Sibling priority ties detected ({} violation{}):\n{}",
                violations.len(),
                if violations.len() == 1 { "" } else { "s" },
                violations.join("\n")
            ));
        }

        self.priority = priorities;
        Ok(())
    }

    /// Returns whether priorities have been set.
    pub fn has_priorities(&self) -> bool {
        !self.priority.is_empty()
    }

    /// Returns the lowest common ancestor of nodes `a` and `b`.
    pub fn lca(&self, a: usize, b: usize) -> usize {
        self.lca_support.lca(a, b)
    }

    /// Returns the lowest common ancestor of nodes `a` and `b`,
    /// treating `None` as the identity element.
    pub fn lca_options(&self, a: Option<usize>, b: Option<usize>) -> Option<usize> {
        if a.is_none() { return b }
        if b.is_none() { return a }
        let (a,b) = (a.unwrap(), b.unwrap());
        Some(self.lca(a,b))
    }

    /// Priority-aware LCA. See the struct-level documentation for a full
    /// description of the priority system.
    ///
    /// All operations are O(1) after construction.
    #[inline]
    pub fn plca(&self, a: usize, b: usize) -> usize {
        if a == b {
            return a;
        }

        let l = self.lca_support.lca(a, b);

        // If no priorities are set, or one node is an ancestor of the
        // other (l == a or l == b), the standard LCA is the answer.
        if self.priority.is_empty() || l == a || l == b {
            return l;
        }

        // Find the children of L on the paths to a and b.
        let depth_l = self.lca_support.depth(l);
        let ca = self.ancestors_by_depth[a][self.lca_support.depth(a) - depth_l - 1];
        let cb = self.ancestors_by_depth[b][self.lca_support.depth(b) - depth_l - 1];

        let pa = self.priority[ca];
        let pb = self.priority[cb];
        if pa < pb {
            a  // ca's subtree wins; a survives
        } else if pb < pa {
            b  // cb's subtree wins; b survives
        } else {
            l  // all siblings share the same priority; standard LCA behavior
        }
    }

    /// Priority-aware LCA with `Option` wrappers, treating `None` as identity.
    /// This is the drop-in replacement for `lca_options` in the fold pattern.
    pub fn plca_options(&self, a: Option<usize>, b: Option<usize>) -> Option<usize> {
        if a.is_none() { return b }
        if b.is_none() { return a }
        Some(self.plca(a.unwrap(), b.unwrap()))
    }

    pub fn n_nodes(&self) -> usize {
        self.n
    }

    pub fn root(&self) -> usize {
        self.root
    }

    /// Returns the parent of `node`. For the root, returns `root` (itself).
    pub fn parent(&self, node: usize) -> usize {
        assert!(node < self.n);
        self.parent[node]
    }

    // --- Serialization ---

    /// Write the tree to `w` in a simple binary format:
    /// each `Vec<usize>` is stored as a little-endian u64 length followed by its elements
    /// as little-endian u64 values.
    pub fn serialize<W: Write>(&self, w: &mut W) -> io::Result<()> {
        write_u64(w, self.n as u64)?;
        write_u64(w, self.root as u64)?;
        write_usize_slice(w, &self.parent)?;
        self.lca_support.serialize(w)?;
        // Priorities: write length then data (length 0 means no priorities).
        write_usize_slice(w, &self.priority)?;
        // Ancestors by depth: write as n vectors.
        write_u64(w, self.ancestors_by_depth.len() as u64)?;
        for anc in &self.ancestors_by_depth {
            write_usize_slice(w, anc)?;
        }
        Ok(())
    }

    /// Load a tree previously written with [`serialize`].
    pub fn load<R: Read>(r: &mut R) -> io::Result<Self> {
        let n = read_u64(r)? as usize;
        let root = read_u64(r)? as usize;
        let parent = read_usize_vec(r)?;
        let lca_support = LcaSupport::load(r)?;
        let priority = read_usize_vec(r)?;
        let n_anc = read_u64(r)? as usize;
        let mut ancestors_by_depth = Vec::with_capacity(n_anc);
        for _ in 0..n_anc {
            ancestors_by_depth.push(read_usize_vec(r)?);
        }
        Ok(LcaTree { n, root, parent, lca_support, priority, ancestors_by_depth })
    }
}

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

    fn binary_tree() -> LcaTree {
        // Complete binary tree, leaves first then internal nodes:
        //         6        <- root (internal)
        //        / \
        //       4   5      <- internal
        //      / \ / \
        //     0  1 2  3   <- leaves
        LcaTree::new(7, vec![(0,4),(1,4),(2,5),(3,5),(4,6),(5,6)]).unwrap()
    }

    fn path_tree(n: usize) -> LcaTree {
        // One leaf (node 0), internal nodes 1..n-1, root = n-1.
        // Chain: 0 -> 1 -> 2 -> ... -> n-1
        let edges = (0..n-1).map(|i| (i, i+1)).collect();
        LcaTree::new(n, edges).unwrap()
    }

    #[test]
    fn single_node() {
        let t = LcaTree::new(1, vec![]).unwrap();
        assert_eq!(t.root(), 0);
        assert_eq!(t.parent(0), 0);
        assert_eq!(t.lca(0, 0), 0);
    }

    #[test]
    fn star_lca() {
        // Leaves: 0,1,2,3  Internal/root: 4
        let t = LcaTree::new(5, vec![(0,4),(1,4),(2,4),(3,4)]).unwrap();
        assert_eq!(t.root(), 4);
        assert_eq!(t.lca(0, 1), 4);
        assert_eq!(t.lca(2, 3), 4);
        assert_eq!(t.lca(0, 4), 4);
        assert_eq!(t.lca(1, 1), 1);
    }

    #[test]
    fn binary_tree_lca() {
        let t = binary_tree();
        assert_eq!(t.root(), 6);

        assert_eq!(t.lca(0, 1), 4);
        assert_eq!(t.lca(2, 3), 5);
        assert_eq!(t.lca(0, 2), 6);
        assert_eq!(t.lca(1, 3), 6);
        assert_eq!(t.lca(0, 4), 4);
        assert_eq!(t.lca(4, 5), 6);
        assert_eq!(t.lca(6, 0), 6);
        assert_eq!(t.lca(3, 3), 3);

        assert_eq!(t.parent(0), 4);
        assert_eq!(t.parent(4), 6);
        assert_eq!(t.parent(6), 6); // root is its own parent
    }

    #[test]
    fn path_tree_lca() {
        // Chain: 0 -> 1 -> 2 -> 3 -> 4 (root)
        let t = path_tree(5);
        assert_eq!(t.root(), 4);

        // LCA is always the shallower (higher-ID) node since it's a chain
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
        let t2 = LcaTree::load(&mut buf.as_slice()).unwrap();

        assert_eq!(t2.n_nodes(), 7);
        assert_eq!(t2.root(), t.root());
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(t2.lca(i, j), t.lca(i, j));
            }
            assert_eq!(t2.parent(i), t.parent(i));
        }
    }

    #[test]
    fn serialize_roundtrip_path() {
        let t = path_tree(10);
        let mut buf = Vec::new();
        t.serialize(&mut buf).unwrap();
        let t2 = LcaTree::load(&mut buf.as_slice()).unwrap();

        assert_eq!(t2.n_nodes(), 10);
        assert_eq!(t2.root(), t.root());
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(t2.lca(i, j), t.lca(i, j));
            }
        }
    }

    // --- Validation error cases ---

    #[test]
    fn error_wrong_edge_count() {
        let err = LcaTree::new(4, vec![(1, 3), (2, 3)]).unwrap_err();
        assert!(err.contains("edges"), "{err}");
    }

    #[test]
    fn error_self_loop() {
        let err = LcaTree::new(3, vec![(0, 2), (2, 2)]).unwrap_err();
        assert!(err.contains("Self-loop"), "{err}");
    }

    #[test]
    fn error_multiple_parents() {
        // Node 0 gets two parents
        let err = LcaTree::new(4, vec![(0, 2), (0, 3), (2, 3)]).unwrap_err();
        assert!(err.contains("multiple parents"), "{err}");
    }

    #[test]
    fn error_disconnected() {
        // Nodes 1,2,3 form a cycle; node 0 is isolated (becomes "root" but unreachable from it)
        let err = LcaTree::new(4, vec![(1, 2), (2, 3), (3, 1)]).unwrap_err();
        assert!(err.contains("disconnected"), "{err}");
    }

    #[test]
    fn error_wrong_leaf_internal_order() {
        // Tree is valid structurally, but internal node 0 has a smaller ID than leaf 1.
        //   0 (internal, root)
        //  / \
        // 1   2  <- leaves, but IDs 1,2 are >= internal ID 0: violation
        let err = LcaTree::new(3, vec![(1, 0), (2, 0)]).unwrap_err();
        assert!(err.contains("leaf") || err.contains("internal"), "{err}");
    }

    #[test]
    fn error_out_of_range_node() {
        let err = LcaTree::new(3, vec![(0, 2), (5, 2)]).unwrap_err();
        assert!(err.contains("out-of-range"), "{err}");
    }

    #[test]
    fn error_zero_nodes() {
        let err = LcaTree::new(0, vec![]).unwrap_err();
        assert!(err.contains("at least one"), "{err}");
    }

    // --- Priority-aware LCA (plca) tests ---

    /// Build a chromosome-inspired hierarchy where leaf siblings share the same
    /// priority (all-same is valid — falls back to standard LCA within each group).
    ///
    ///           6 (root)
    ///          / \
    ///   4 (arm)   5 (centromere)    <- centromere has higher priority
    ///    / \       / \
    ///   0   1     2   3
    ///  (p) (q)  (cL) (cR)
    fn chromosome_tree_with_ties() -> LcaTree {
        let mut t = LcaTree::new(7, vec![
            (0, 4), (1, 4),
            (2, 5), (3, 5),
            (4, 6), (5, 6),
        ]).unwrap();
        // Leaves under arm all share priority 1 (valid: all same).
        // Leaves under centromere all share priority 1 (valid: all same).
        // arm=2, centromere=1 under root (valid: all distinct).
        //                 0  1  2  3  4(arm) 5(cen) 6(root)
        t.set_priorities(vec![1, 1, 1, 1, 2,     1,     99]).unwrap();
        t
    }

    #[test]
    fn plca_without_priorities_equals_lca() {
        let t = binary_tree(); // no priorities set
        assert!(!t.has_priorities());
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(t.plca(i, j), t.lca(i, j),
                    "plca should equal lca when no priorities are set");
            }
        }
    }

    #[test]
    fn plca_centromere_wins_over_arm() {
        let t = chromosome_tree_with_ties();
        // p_arm (0) vs cen_left (2): LCA would be root (6), but centromere (pri 1)
        // beats arm (pri 2) → cen_left survives
        assert_eq!(t.plca(0, 2), 2);
        assert_eq!(t.plca(1, 3), 3);
        assert_eq!(t.plca(0, 3), 3);
        assert_eq!(t.plca(1, 2), 2);
    }

    #[test]
    fn plca_tie_returns_lca() {
        let t = chromosome_tree_with_ties();
        // p_arm vs q_arm: same priority under arm → tie → returns arm
        assert_eq!(t.plca(0, 1), 4);
        // cen_left vs cen_right: same priority under centromere → tie → returns centromere
        assert_eq!(t.plca(2, 3), 5);
    }

    /// Same chromosome tree but with distinct sub-priorities to avoid ties.
    fn chromosome_tree_no_ties() -> LcaTree {
        let mut t = LcaTree::new(7, vec![
            (0, 4), (1, 4),
            (2, 5), (3, 5),
            (4, 6), (5, 6),
        ]).unwrap();
        //                 0(p) 1(q) 2(cL) 3(cR) 4(arm) 5(cen) 6(root)
        t.set_priorities(vec![1,   2,   1,    2,    2,     1,     99]).unwrap();
        t
    }

    #[test]
    fn plca_no_ties_centromere_wins() {
        let t = chromosome_tree_no_ties();
        // p_arm vs cen_left: centromere subtree wins → cen_left
        assert_eq!(t.plca(0, 2), 2);
        // q_arm vs cen_right: centromere subtree wins → cen_right
        assert_eq!(t.plca(1, 3), 3);
    }

    #[test]
    fn plca_no_ties_within_subtree() {
        let t = chromosome_tree_no_ties();
        // p_arm (pri 1) vs q_arm (pri 2) under arm: p_arm wins
        assert_eq!(t.plca(0, 1), 0);
        // cen_left (pri 1) vs cen_right (pri 2) under centromere: cen_left wins
        assert_eq!(t.plca(2, 3), 2);
    }

    #[test]
    fn plca_is_symmetric() {
        let t = chromosome_tree_no_ties();
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(t.plca(i, j), t.plca(j, i),
                    "PLCA({i},{j}) != PLCA({j},{i})");
            }
        }
    }

    #[test]
    fn plca_with_self_is_self() {
        let t = chromosome_tree_no_ties();
        for i in 0..7 {
            assert_eq!(t.plca(i, i), i);
        }
    }

    #[test]
    fn plca_associativity_no_ties() {
        let t = chromosome_tree_no_ties();
        for a in 0..7 {
            for b in 0..7 {
                for c in 0..7 {
                    let left = t.plca(t.plca(a, b), c);
                    let right = t.plca(a, t.plca(b, c));
                    assert_eq!(left, right,
                        "Associativity failed: plca(plca({a},{b}),{c})={left} != plca({a},plca({b},{c}))={right}");
                }
            }
        }
    }

    #[test]
    fn plca_ancestor_subsumes_descendant() {
        let t = chromosome_tree_no_ties();
        // arm (4) vs p_arm (0): arm is ancestor → returns arm
        assert_eq!(t.plca(4, 0), 4);
        assert_eq!(t.plca(0, 4), 4);
        // root (6) vs anything → root
        for i in 0..7 {
            assert_eq!(t.plca(6, i), 6);
        }
    }

    #[test]
    fn plca_fold_order_independent() {
        let t = chromosome_tree_no_ties();
        // Simulate folding a sequence of colors in different orders.
        // p_arm, q_arm, cen_left → should always give cen_left (centromere wins)
        let colors = vec![0, 1, 2]; // p_arm, q_arm, cen_left
        let forward = colors.iter().copied().reduce(|a, b| t.plca(a, b)).unwrap();
        let reverse = colors.iter().rev().copied().reduce(|a, b| t.plca(a, b)).unwrap();
        assert_eq!(forward, reverse, "Fold should be order-independent");
        assert_eq!(forward, 2, "cen_left should win");
    }

    #[test]
    fn plca_options_identity() {
        let t = chromosome_tree_no_ties();
        // None is identity
        assert_eq!(t.plca_options(None, Some(3)), Some(3));
        assert_eq!(t.plca_options(Some(3), None), Some(3));
        assert_eq!(t.plca_options(None, None), None);
        // Normal case
        assert_eq!(t.plca_options(Some(0), Some(2)), Some(2));
    }

    #[test]
    fn set_priorities_accepts_all_same() {
        let mut t = binary_tree();
        //         6
        //        / \
        //       4   5      <- same priority: all-same is valid
        //      / \ / \
        //     0  1 2  3    <- same priority under each parent: all-same is valid
        t.set_priorities(vec![1, 1, 1, 1, 1, 1, 99]).unwrap();
        // plca should behave like lca when siblings share priorities
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(t.plca(i, j), t.lca(i, j),
                    "plca should equal lca when all siblings share priority");
            }
        }
    }

    #[test]
    fn set_priorities_rejects_mixed() {
        let mut t = binary_tree();
        //         6
        //        / \
        //       4   5      <- distinct priorities (1, 2): valid
        //      / \ / \
        //     0  1 2  3
        // Give 0 and 1 priorities 1 and 1 (all-same: valid),
        // but give 2 and 3 priorities 1 and 2 (all-distinct: valid).
        // This should pass.
        t.set_priorities(vec![1, 1, 1, 2, 1, 2, 99]).unwrap();

        // Now try a mixed case: 3 siblings under one parent, two share a
        // priority and one differs.
        // Star: 0,1,2,3 are leaves under root 4.
        let mut star = LcaTree::new(5, vec![(0,4),(1,4),(2,4),(3,4)]).unwrap();
        //                          0  1  2  3  4(root)
        let err = star.set_priorities(vec![1, 1, 2, 3, 99]).unwrap_err();
        assert!(err.contains("duplicate priority 1"), "{err}");
        assert!(err.contains("0") && err.contains("1"), "{err}");
    }

    #[test]
    fn set_priorities_reports_all_violations() {
        // Binary tree where BOTH sibling groups have mixed priorities.
        //         6
        //        / \
        //       4   5
        //      / \ / \
        //     0  1 2  3
        // Under node 4: give 0 and 1 priorities 1 and 1 — all-same, valid.
        // Under node 5: give 2 and 3 priorities 1 and 1 — all-same, valid.
        // Under node 6: give 4 and 5 priorities 1 and 1 — all-same, valid.
        // That's fine. Now make it mixed:

        // Star with 5 children, some sharing priorities, some not.
        let mut star = LcaTree::new(6, vec![(0,5),(1,5),(2,5),(3,5),(4,5)]).unwrap();
        //                          0  1  2  3  4  5(root)
        let err = star.set_priorities(vec![1, 1, 2, 3, 3, 99]).unwrap_err();
        assert!(err.contains("2 violations"), "{err}");
        // Should report both the [0,1] group and the [3,4] group
        assert!(err.contains("duplicate priority 1"), "{err}");
        assert!(err.contains("duplicate priority 3"), "{err}");
    }

    #[test]
    fn set_priorities_wrong_length() {
        let mut t = binary_tree();
        let err = t.set_priorities(vec![1, 2, 3]).unwrap_err();
        assert!(err.contains("Expected 7"), "{err}");
    }

    #[test]
    fn serialize_roundtrip_with_priorities() {
        let t = chromosome_tree_no_ties();
        let mut buf = Vec::new();
        t.serialize(&mut buf).unwrap();
        let t2 = LcaTree::load(&mut buf.as_slice()).unwrap();

        assert_eq!(t2.n_nodes(), 7);
        assert_eq!(t2.root(), t.root());
        assert!(t2.has_priorities());
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(t2.lca(i, j), t.lca(i, j));
                assert_eq!(t2.plca(i, j), t.plca(i, j),
                    "PLCA({i},{j}) mismatch after roundtrip");
            }
            assert_eq!(t2.parent(i), t.parent(i));
        }
    }
}
