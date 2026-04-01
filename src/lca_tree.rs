use std::io::{self, Read, Write};
use crate::lca_support::LcaSupport;

/// A rooted tree where leaves occupy IDs 0..num_leaves and internal nodes occupy
/// IDs num_leaves..n. The root is an internal node (or node 0 for a single-node tree).
/// `parent[root] == root` (root is its own parent).
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct LcaTree {
    n: usize,
    root: usize,
    /// parent[i] is the parent of node i; parent[root] == root.
    parent: Vec<usize>,
    lca_support: LcaSupport,
}

impl LcaTree {
    /// Construct a tree from `n` nodes and `edges`, where each edge `(child, parent)`
    /// points toward the root. Returns an error if these conditions are not met or
    /// if the edges do not form a valid tree.
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

        let lca_support = LcaSupport::new(n, edges).unwrap();

        Ok(LcaTree { n, root, parent, lca_support })
    }

    /// Returns the lowest common ancestor of nodes `a` and `b`.
    pub fn lca(&self, a: usize, b: usize) -> usize {
        self.lca_support.lca(a, b)
    }

    /// Returns true if `candidate` is an ancestor-or-equal of `node`.
    /// Uses LCA for O(1) query: candidate is an ancestor of node iff lca(node, candidate) == candidate.
    pub fn is_ancestor(&self, node: usize, candidate: usize) -> bool {
        self.lca(node, candidate) == candidate
    }

    /// Returns an iterator over all strict ancestors of `node` (not including `node` itself),
    /// from parent up to the root, by walking the parent array.
    pub fn ancestors(&self, node: usize) -> impl Iterator<Item = usize> + '_ {
        let mut cur = Some(node);
        std::iter::from_fn(move || {
            let val = cur?;
            cur = if self.parent[val] == val { None } else { Some(self.parent[val]) };
            Some(val)
        })
    }

    /// Returns the lowest common ancestor of nodes `a` and `b`.
    pub fn lca_options(&self, a: Option<usize>, b: Option<usize>) -> Option<usize> {
        if a.is_none() { return b }
        if b.is_none() { return a }
        let (a,b) = (a.unwrap(), b.unwrap());
        Some(self.lca(a,b))
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
        Ok(())
    }

    /// Load a tree previously written with [`serialize`].
    pub fn load<R: Read>(r: &mut R) -> io::Result<Self> {
        let n = read_u64(r)? as usize;
        let root = read_u64(r)? as usize;
        let parent = read_usize_vec(r)?;
        let lca_support = LcaSupport::load(r)?;
        Ok(LcaTree { n, root, parent, lca_support })
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
    fn error_out_of_range_node() {
        let err = LcaTree::new(3, vec![(0, 2), (5, 2)]).unwrap_err();
        assert!(err.contains("out-of-range"), "{err}");
    }

    #[test]
    fn error_zero_nodes() {
        let err = LcaTree::new(0, vec![]).unwrap_err();
        assert!(err.contains("at least one"), "{err}");
    }

    // --- is_ancestor ---

    #[test]
    fn is_ancestor_binary_tree() {
        let t = binary_tree();
        // root (6) is ancestor of everyone
        for i in 0..7 {
            assert!(t.is_ancestor(i, 6), "6 should be ancestor of {i}");
        }
        // internal node 4 is ancestor of its leaves 0,1 and itself, not of 2,3,5
        assert!(t.is_ancestor(0, 4));
        assert!(t.is_ancestor(1, 4));
        assert!(t.is_ancestor(4, 4)); // reflexive
        assert!(!t.is_ancestor(2, 4));
        assert!(!t.is_ancestor(3, 4));
        assert!(!t.is_ancestor(5, 4));
        // leaves are not ancestors of anything except themselves
        assert!(t.is_ancestor(0, 0));
        assert!(!t.is_ancestor(1, 0));
        assert!(!t.is_ancestor(6, 0));
    }

    #[test]
    fn is_ancestor_path_tree() {
        // Chain: 0 -> 1 -> 2 -> 3 -> 4 (root)
        let t = path_tree(5);
        // Every node is an ancestor of 0
        for anc in 0..5 {
            assert!(t.is_ancestor(0, anc), "{anc} should be ancestor of 0");
        }
        // 0 is only an ancestor of itself
        for desc in 1..5 {
            assert!(!t.is_ancestor(desc, 0), "0 should not be ancestor of {desc}");
        }
        // 2 is an ancestor of 1 and 0, not of 3 or 4's subtree peers
        assert!(t.is_ancestor(1, 2));
        assert!(t.is_ancestor(0, 2));
        assert!(!t.is_ancestor(3, 2));
        assert!(!t.is_ancestor(4, 2));
    }

    // --- ancestors ---

    #[test]
    fn ancestors_binary_tree() {
        let t = binary_tree();
        // Leaf 0: self + ancestors should be 0, 4, 6
        let ancs: Vec<usize> = t.ancestors(0).collect();
        assert_eq!(ancs, vec![0, 4, 6]);
        // Leaf 2: self + ancestors should be 2, 5, 6
        let ancs: Vec<usize> = t.ancestors(2).collect();
        assert_eq!(ancs, vec![2, 5, 6]);
        // Internal node 4: self + ancestors should be 4, 6
        let ancs: Vec<usize> = t.ancestors(4).collect();
        assert_eq!(ancs, vec![4, 6]);
        // Root 6: just itself
        let ancs: Vec<usize> = t.ancestors(6).collect();
        assert_eq!(ancs, vec![6]);
    }

    #[test]
    fn ancestors_path_tree() {
        // Chain: 0 -> 1 -> 2 -> 3 -> 4 (root)
        let t = path_tree(5);
        let ancs: Vec<usize> = t.ancestors(0).collect();
        assert_eq!(ancs, vec![0, 1, 2, 3, 4]);
        let ancs: Vec<usize> = t.ancestors(3).collect();
        assert_eq!(ancs, vec![3, 4]);
        let ancs: Vec<usize> = t.ancestors(4).collect();
        assert_eq!(ancs, vec![4]);
    }
}
