use std::collections::HashMap;
use crate::lca_tree::LcaTree;

/// Construction-time priority-aware LCA. Borrows an [`LcaTree`] and augments
/// it with per-node priorities plus an ancestor table needed for O(1) `plca`.
/// Not serialized: build it, fold with it, drop it before writing the index.
///
/// With standard LCA using the chromosome hierarchy, a k-mer found in both
/// "chr1" and "chrM" would be labeled as their common ancestor (e.g. "root"),
/// losing specificity. Using priorities, the user can declare that
/// "mitochondrial" is more interesting than "autosomal" for their use case.
///
/// Each node is assigned an integer priority (lower value = higher priority).
/// When computing `plca(a, b)`:
/// - If `a == b` or one is an ancestor of the other, returns the standard LCA.
/// - Otherwise, let `ca` and `cb` be the children of `LCA(a, b)` that are
///   ancestors of `a` and `b` respectively. The side with the smaller
///   priority value wins (`a` or `b` is returned); on tie, the standard LCA
///   is returned.
///
/// For any group of siblings, priorities must be either all equal or all
/// distinct. This ensures `plca` is associative and safe for incremental
/// folds. Mixed groups are rejected by [`PriorityLca::new`].
#[derive(Debug)]
pub struct PriorityLca<'a> {
    tree: &'a LcaTree,
    /// Lower value = higher priority.
    priority: Vec<usize>,
    /// `ancestors_by_depth[v]` = `[v, parent(v), grandparent(v), ..., root]`.
    /// Used to answer "child of ancestor A on path to descendant D" in O(1):
    ///   `ancestors_by_depth[D][depth(D) - depth(A) - 1]`.
    ancestors_by_depth: Vec<Vec<usize>>,
}

impl<'a> PriorityLca<'a> {
    /// Build a priority-aware LCA over `tree` with one priority per node.
    /// Returns an error if `priority.len() != tree.n_nodes()` or any sibling
    /// group has a mix of shared and distinct priorities.
    pub fn new(tree: &'a LcaTree, priority: Vec<usize>) -> Result<Self, String> {
        let n = tree.n_nodes();
        if priority.len() != n {
            return Err(format!(
                "Expected {} priorities, got {}",
                n,
                priority.len()
            ));
        }

        // Reconstruct children lists from parent pointers.
        let root = tree.root();
        let mut children = vec![vec![]; n];
        for v in 0..n {
            if v != root {
                children[tree.parent(v)].push(v);
            }
        }

        validate_sibling_priorities(&children, &priority)?;

        // ancestors_by_depth[v] = [v, parent(v), ..., root]
        let mut ancestors_by_depth = vec![vec![]; n];
        for v in 0..n {
            let mut anc = Vec::with_capacity(tree.depth(v) + 1);
            let mut cur = v;
            loop {
                anc.push(cur);
                if cur == root { break; }
                cur = tree.parent(cur);
            }
            ancestors_by_depth[v] = anc;
        }

        Ok(Self { tree, priority, ancestors_by_depth })
    }

    /// Priority-aware LCA. O(1) after construction.
    #[inline]
    pub fn plca(&self, a: usize, b: usize) -> usize {
        if a == b {
            return a;
        }

        let l = self.tree.lca(a, b);
        if l == a || l == b {
            return l;
        }

        let dl = self.tree.depth(l);
        let ca = self.ancestors_by_depth[a][self.tree.depth(a) - dl - 1];
        let cb = self.ancestors_by_depth[b][self.tree.depth(b) - dl - 1];

        use std::cmp::Ordering::*;
        match self.priority[ca].cmp(&self.priority[cb]) {
            Less    => a,
            Greater => b,
            Equal   => l,
        }
    }

    /// Priority-aware LCA with `Option` wrappers, treating `None` as the
    /// identity element. Drop-in replacement for `LcaTree::lca_options` in
    /// the fold pattern.
    pub fn plca_options(&self, a: Option<usize>, b: Option<usize>) -> Option<usize> {
        match (a, b) {
            (None, x) | (x, None) => x,
            (Some(a), Some(b))    => Some(self.plca(a, b)),
        }
    }

    pub fn tree(&self) -> &LcaTree {
        self.tree
    }
}

/// Validate that every sibling group's priorities are either all equal or
/// all distinct. Collects and reports every violation.
fn validate_sibling_priorities(
    children: &[Vec<usize>],
    priority: &[usize],
) -> Result<(), String> {
    let mut violations = Vec::new();

    for (v, kids) in children.iter().enumerate() {
        if kids.len() < 2 {
            continue;
        }
        let first = priority[kids[0]];
        if kids.iter().all(|&c| priority[c] == first) {
            continue; // all-same: valid
        }

        // Not all-same: any duplicates are violations.
        let mut by_priority: HashMap<usize, Vec<usize>> = HashMap::new();
        for &c in kids {
            by_priority.entry(priority[c]).or_default().push(c);
        }
        let mut dup_priorities: Vec<usize> = by_priority
            .iter()
            .filter(|(_, ns)| ns.len() > 1)
            .map(|(p, _)| *p)
            .collect();
        dup_priorities.sort();
        for pri in dup_priorities {
            let nodes = &by_priority[&pri];
            let node_list: Vec<String> = nodes.iter().map(|n| n.to_string()).collect();
            violations.push(format!(
                "  Children of node {v} with duplicate priority {pri}: [{}]. \
                 Either assign them distinct priorities or add a named grouping \
                 node as their parent at priority {pri}.",
                node_list.join(", ")
            ));
        }
    }

    if violations.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "Sibling priority ties detected ({} violation{}):\n{}",
            violations.len(),
            if violations.len() == 1 { "" } else { "s" },
            violations.join("\n")
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn binary_tree() -> LcaTree {
        //         6        <- root
        //        / \
        //       4   5
        //      / \ / \
        //     0  1 2  3
        LcaTree::new(7, vec![(0,4),(1,4),(2,5),(3,5),(4,6),(5,6)]).unwrap()
    }

    /// Chromosome-inspired hierarchy with distinct sub-priorities.
    ///
    ///           6 (root)
    ///          / \
    ///   4 (arm)   5 (centromere)    <- centromere has higher priority (1 < 2)
    ///    / \       / \
    ///   0   1     2   3
    fn chromosome_tree_no_ties<'a>(tree: &'a LcaTree) -> PriorityLca<'a> {
        //                 0(p) 1(q) 2(cL) 3(cR) 4(arm) 5(cen) 6(root)
        PriorityLca::new(tree, vec![1,   2,   1,    2,    2,     1,     99]).unwrap()
    }

    fn chromosome_tree_with_ties<'a>(tree: &'a LcaTree) -> PriorityLca<'a> {
        //                 0  1  2  3  4(arm) 5(cen) 6(root)
        PriorityLca::new(tree, vec![1, 1, 1, 1, 2,     1,     99]).unwrap()
    }

    #[test]
    fn plca_centromere_wins_over_arm() {
        let t = binary_tree();
        let p = chromosome_tree_with_ties(&t);
        assert_eq!(p.plca(0, 2), 2);
        assert_eq!(p.plca(1, 3), 3);
        assert_eq!(p.plca(0, 3), 3);
        assert_eq!(p.plca(1, 2), 2);
    }

    #[test]
    fn plca_tie_returns_lca() {
        let t = binary_tree();
        let p = chromosome_tree_with_ties(&t);
        assert_eq!(p.plca(0, 1), 4);
        assert_eq!(p.plca(2, 3), 5);
    }

    #[test]
    fn plca_no_ties_centromere_wins() {
        let t = binary_tree();
        let p = chromosome_tree_no_ties(&t);
        assert_eq!(p.plca(0, 2), 2);
        assert_eq!(p.plca(1, 3), 3);
    }

    #[test]
    fn plca_no_ties_within_subtree() {
        let t = binary_tree();
        let p = chromosome_tree_no_ties(&t);
        assert_eq!(p.plca(0, 1), 0);
        assert_eq!(p.plca(2, 3), 2);
    }

    #[test]
    fn plca_is_symmetric() {
        let t = binary_tree();
        let p = chromosome_tree_no_ties(&t);
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(p.plca(i, j), p.plca(j, i));
            }
        }
    }

    #[test]
    fn plca_with_self_is_self() {
        let t = binary_tree();
        let p = chromosome_tree_no_ties(&t);
        for i in 0..7 {
            assert_eq!(p.plca(i, i), i);
        }
    }

    #[test]
    fn plca_associativity_no_ties() {
        let t = binary_tree();
        let p = chromosome_tree_no_ties(&t);
        for a in 0..7 {
            for b in 0..7 {
                for c in 0..7 {
                    assert_eq!(p.plca(p.plca(a, b), c), p.plca(a, p.plca(b, c)));
                }
            }
        }
    }

    #[test]
    fn plca_ancestor_subsumes_descendant() {
        let t = binary_tree();
        let p = chromosome_tree_no_ties(&t);
        assert_eq!(p.plca(4, 0), 4);
        assert_eq!(p.plca(0, 4), 4);
        for i in 0..7 {
            assert_eq!(p.plca(6, i), 6);
        }
    }

    #[test]
    fn plca_fold_order_independent() {
        let t = binary_tree();
        let p = chromosome_tree_no_ties(&t);
        let colors = vec![0, 1, 2];
        let fwd = colors.iter().copied().reduce(|a, b| p.plca(a, b)).unwrap();
        let rev = colors.iter().rev().copied().reduce(|a, b| p.plca(a, b)).unwrap();
        assert_eq!(fwd, rev);
        assert_eq!(fwd, 2);
    }

    #[test]
    fn plca_options_identity() {
        let t = binary_tree();
        let p = chromosome_tree_no_ties(&t);
        assert_eq!(p.plca_options(None, Some(3)), Some(3));
        assert_eq!(p.plca_options(Some(3), None), Some(3));
        assert_eq!(p.plca_options(None, None), None);
        assert_eq!(p.plca_options(Some(0), Some(2)), Some(2));
    }

    #[test]
    fn new_accepts_all_same() {
        let t = binary_tree();
        let p = PriorityLca::new(&t, vec![1, 1, 1, 1, 1, 1, 99]).unwrap();
        // plca should equal lca when all sibling groups share priority
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(p.plca(i, j), t.lca(i, j));
            }
        }
    }

    #[test]
    fn new_rejects_mixed() {
        let star = LcaTree::new(5, vec![(0,4),(1,4),(2,4),(3,4)]).unwrap();
        let err = PriorityLca::new(&star, vec![1, 1, 2, 3, 99]).unwrap_err();
        assert!(err.contains("duplicate priority 1"), "{err}");
    }

    #[test]
    fn new_reports_all_violations() {
        let star = LcaTree::new(6, vec![(0,5),(1,5),(2,5),(3,5),(4,5)]).unwrap();
        let err = PriorityLca::new(&star, vec![1, 1, 2, 3, 3, 99]).unwrap_err();
        assert!(err.contains("2 violations"), "{err}");
        assert!(err.contains("duplicate priority 1"), "{err}");
        assert!(err.contains("duplicate priority 3"), "{err}");
    }

    #[test]
    fn new_wrong_length() {
        let t = binary_tree();
        let err = PriorityLca::new(&t, vec![1, 2, 3]).unwrap_err();
        assert!(err.contains("Expected 7"), "{err}");
    }
}
