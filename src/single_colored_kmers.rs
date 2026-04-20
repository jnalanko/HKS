use std::io::{Read, Write};
use std::ops::Range;

use bitvec::prelude::*;
use bitvec_sds::traits::RandomAccessU32;
//use bitvec_sds::wavelet_tree::{SelectSupportBoth, WaveletTree};
use sbwt::{ContractLeft, LcsArray, MatchingStatisticsIterator, SbwtIndex, StreamingIndex, SubsetMatrix};
use crate::color_storage::SimpleColorStorage;
use crate::lca_tree::LcaTree;
use crate::traits::*;

#[derive(Debug, Clone)]
pub struct HksIndex<L: ContractLeft + Clone + MySerialize + From<LcsArray>, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: L,
    feature_set: Labeling<C>,
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> HksIndex<L, C> {

    pub fn feature_set(&self) -> &Labeling<C> {
        &self.feature_set
    }

    pub fn into_parts(self) -> (SbwtIndex<SubsetMatrix>, L, Labeling<C>) {
        (self.sbwt, self.lcs, self.feature_set)
    }

    pub fn rename_labels(&mut self, new_names: Vec<String>) {
        self.feature_set.hierarchy.rename_labels(new_names);
    }

    pub fn k(&self) -> usize {
        self.sbwt.k()
    }

    pub fn sbwt(&self) -> &SbwtIndex<SubsetMatrix> {
        &self.sbwt
    }

    pub fn lcs(&self) -> &L {
        &self.lcs
    }

    fn base_serialization_magic() -> [u8; 4] {
        [17, 42, 191, 203]
    }

    fn base_serialization_version() -> u32 {
        7_u32
    }

    pub fn serialize_base(&self, mut out: &mut impl Write) {
        out.write_all(&Self::base_serialization_magic()).unwrap();
        out.write_all(&Self::base_serialization_version().to_le_bytes()).unwrap();
        self.sbwt.serialize(out).unwrap();
        self.lcs.serialize(out);
    }

    /// Load sbwt and lcs from a base index file. Use `from_parts` to combine with a loaded
    /// `FeatureSet` into a full `HksIndex`.
    pub fn load_base(mut input: &mut impl Read) -> (SbwtIndex<SubsetMatrix>, L) {
        let mut magic = [0_u8; 4];
        input.read_exact(&mut magic).unwrap();
        if magic != Self::base_serialization_magic() {
            panic!("Error loading index: invalid file format (magic constant mismatch)");
        }

        let mut version_bytes = [0_u8; 4];
        input.read_exact(&mut version_bytes).unwrap();
        let version = u32::from_le_bytes(version_bytes);
        if version != Self::base_serialization_version() {
            panic!("Error loading index: wrong file format version number (found {}, expected {})", version, Self::base_serialization_version());
        }

        let sbwt = SbwtIndex::<sbwt::SubsetMatrix>::load(input).unwrap();
        let lcs = *L::load(input);
        (sbwt, lcs)
    }

    pub fn n_kmers(&self) -> usize {
        self.sbwt.n_kmers()
    }

    fn color_run_stats(&self) -> (usize, usize, f64) { // (min, max, mean)
        let n = self.sbwt.n_sets();
        if n == 0 { return (0, 0, 0.0); }

        let mut min = usize::MAX;
        let mut max = 0_usize;
        let mut n_runs = 0_usize;
        let mut prev = self.get_color(0);
        let mut run_len = 1_usize;
        for i in 1..n {
            let cur = self.get_color(i);
            if cur == prev {
                run_len += 1;
            } else {
                min = min.min(run_len);
                max = max.max(run_len);
                n_runs += 1;
                run_len = 1;
                prev = cur;
            }
        }
        min = min.min(run_len);
        max = max.max(run_len);
        n_runs += 1;
        let mean = n as f64 / n_runs as f64;
        (min, max, mean)
    }

    pub fn color_stats(&self) -> ColorStats {
        let mut uncolored = 0_usize;
        let mut colored = 0_usize;
        let mut color_counts = vec![0_usize; self.feature_set.hierarchy.n_nodes()];
        for i in 0..self.sbwt.n_sets() {
            match self.get_color(i) {
                Some(id) => { colored += 1; color_counts[id] += 1; },
                None => uncolored += 1,
            }
        }
        let (color_run_min, color_run_max, color_run_mean) = self.color_run_stats();
        ColorStats { colored, uncolored, color_run_min, color_run_max, color_run_mean, color_counts }
    }

    pub fn get_color(&self, colex: usize) -> Option<usize> {
        assert!(colex < self.sbwt.n_sets());
        self.feature_set.color_assignments.get_color(colex)
    }

    pub fn get_color_of_range(&self, colex_range: Range<usize>) -> Option<usize> {
        assert!(colex_range.end <= self.sbwt.n_sets());
        let fs = &self.feature_set;
        fs.color_assignments.get_color_of_range(colex_range, fs.hierarchy.tree())
    }

    // Returns an iterator giving the color of each of the n-k+1 k-mers of the query.
    // k must be less or equal to the k in the SBWT index.
    // If query is shorter than k, returns an empty iterator.
    pub fn lookup_kmers<'a, 'b>(&'a self, query: &'b [u8], k: usize) -> KmerLookupIterator<'a, 'b, L, C>{
        assert!(k <= self.sbwt.k());
        let si = StreamingIndex {
            extend_right: &self.sbwt,
            contract_left: &self.lcs,
            n: self.sbwt.n_sets(),
            k: self.sbwt.k(),
        };

        //let mut ms_iter = si.bounded_matching_statistics_iter(query, k);
        let mut ms_iter = si.matching_statistics_iter(query);

        // Skip over the first k-1 positions
        for _ in 0..k-1 {
            ms_iter.next(); // If the iterator ends early, will keep returning None
        }
        KmerLookupIterator { matching_stats_iter: ms_iter, index: self, query_pattern_length: k, feature_set: &self.feature_set }
    }

    /// Build a new index from already-converted parts. Use after `new_with_feature_set`
    /// during construction, or after `load_base` + `FeatureSet::load_from_file` during loading.
    pub fn from_parts(sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: L, feature_set: Labeling<C>) -> Self {
        HksIndex::<L, C> { sbwt, lcs, feature_set }
    }

    /// Build a new index from raw construction outputs. Converts the `LcsArray` to `L`.
    pub fn new_with_feature_set(sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: sbwt::LcsArray, feature_set: Labeling<C>) -> Self {
        log::info!("Indexing LCS array");
        let lcs_index = L::from(lcs);
        log::info!("Color structure construction complete");
        HksIndex::<L, C> { sbwt, lcs: lcs_index, feature_set }
    }

    pub fn n_sbwt_sets(&self) -> usize {
        self.sbwt.n_sets()
    }

    pub fn build_sbwt_select(&mut self) {
        self.sbwt.build_select();
    }

    // S is the s-mer length, s <= k
    // Returns a vector of length equal to the number of colors
    // The i-th element in the vector is the number of s-mer assigned
    // to color i. NOTE: build_sbwt_select() must have been called before running this.
    pub fn node_stats(&self, s: usize, dummy_marks: &BitSlice) -> Vec<usize> {
        let feature_set = &self.feature_set;
        let mut counts = vec![0; feature_set.hierarchy.n_nodes()];
        assert!(s <= self.sbwt.k());
        let n = self.n_sbwt_sets();

        // Sweep through every maximal run of positions whose consecutive LCS >= s
        // (i.e. all k-mers in the run share a common s-mer). Compute the LCA of all
        // colors in the run.
        let mut run_start = 0usize;
        for run_end in 1..=n {
            let run_continues = run_end < n && self.lcs.get_lcs(run_end) >= s;
            if !run_continues {
                // Run is run_start..colex
                let mut lca: Option<usize> = None;
                if run_end - run_start == 1 && dummy_marks[run_start] {
                    // We must only count this if the dummy has length at least s.
                    let dummy_len = self.sbwt.access_kmer(run_start).iter().filter(|c| **c != b'$').count();
                    if dummy_len >= s {
                        lca = feature_set.hierarchy.tree().lca_options(lca, feature_set.color_assignments.get_color(run_start));
                    }
                } else {
                    // Since the length of the range is at least 2, all s-mers in the range
                    // are dollar-free: otherwise we would have a duplicate dummy.
                    for pos in run_start..run_end {
                        lca = feature_set.hierarchy.tree().lca_options(lca, feature_set.color_assignments.get_color(pos));
                    }
                }
                if let Some(x) = lca {
                    counts[x] += 1;
                }
                run_start = run_end;
            }
        }

        counts
    }
}

pub struct KmerLookupIterator<'a, 'b, L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    // This iterator should be initialized so that the first k-1 MS values are skipped
    matching_stats_iter: MatchingStatisticsIterator<'a, 'b, SbwtIndex::<SubsetMatrix>, L>,
    index: &'a HksIndex<L, C>,
    query_pattern_length: usize,
    feature_set: &'a Labeling<C>, // Cached to avoid re-fetching on every call
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> Iterator for KmerLookupIterator<'_, '_, L, C> {
    type Item = Option<usize>; // Color id of k-mer

    fn next(&mut self) -> Option<Self::Item> {
        let (len, range) = self.matching_stats_iter.next()?;

        if len >= self.query_pattern_length {
            // k-mer is found in the sbwt
            assert!(range.len() > 0);
            if self.query_pattern_length == self.index.k() {
                assert!(range.len() == 1);
                // No interval expansion or LCA needed, just return the color. Saves some work.
                return Some(self.index.get_color(range.start))
            }

            let lcs = &self.index.lcs;
            let label_tree = &self.feature_set.hierarchy.tree;
            let root_id = label_tree.root();

            let mut color = self.index.get_color_of_range(range.clone());
            if color == Some(root_id) { return Some(color) }

            // Expand the interval as long as it has at most single color and the
            // LCS is at least the query pattern length.
            let (mut new_start, mut new_end) = (range.start, range.end);
            //eprintln!("Expanding from {}..{}, (len {})", new_start, new_end, new_end - new_start);
            while new_start > 0 && lcs.get_lcs(new_start) >= self.query_pattern_length {
                new_start -= 1;
                color = label_tree.lca_options(color, self.index.get_color(new_start));
                if color == Some(root_id) { return Some(color) } // This is a Some(Some(color)). Means that the iterator produced something.
            }
            let n = self.index.sbwt.n_sets();
            while new_end < n && lcs.get_lcs(new_end) >= self.query_pattern_length {
                color = label_tree.lca_options(color, self.index.get_color(new_end));
                if color == Some(root_id) { return Some(color) } // This is a Some(Some(color)). Means that the iterator produced something.
                new_end += 1;
            }

            Some(color) // This is a Some(Some(color)). Means that the iterator produced something.
        } else {
            Some(None) // Iterator not finished but the k-mer is not found -> no color
        }
    }
}

/// A rooted color hierarchy tree together with a name for every node.
/// Leaves occupy IDs 0..n_leaves; internal nodes (including root) occupy n_leaves..n_nodes.
/// `names[i]` is the name of node `i`.
#[derive(Debug, Clone)]
pub struct ColorHierarchy {
    tree: LcaTree,
    names: Vec<String>,
}

impl ColorHierarchy {
    /// Default constructor: creates a star topology with all `leaf_names` as children
    /// of a new root node named "root". No leaf name may be "root".
    pub fn new_star(leaf_names: Vec<String>) -> Self {
        for name in &leaf_names {
            assert!(name != "root", "label name 'root' is reserved");
        }
        let mut names = leaf_names;
        names.push("root".to_string());
        let n = names.len();
        let edges = (0..n - 1).map(|i| (i, n - 1)).collect();
        let tree = LcaTree::new(n, edges).expect("star hierarchy construction cannot fail");
        Self { tree, names }
    }

    /// Constructs a hierarchy from a pre-built `LcaTree` and the names of all nodes.
    pub fn with_tree(tree: LcaTree, node_names: Vec<String>) -> Self {
        assert_eq!(node_names.len(), tree.n_nodes(), "names must have one entry per tree node");
        Self { tree, names: node_names }
    }

    pub fn tree(&self) -> &LcaTree {
        &self.tree
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn n_nodes(&self) -> usize {
        self.tree.n_nodes()
    }

    pub fn root(&self) -> usize {
        self.tree.root()
    }

    pub fn rename_labels(&mut self, new_names: Vec<String>) {
        assert_eq!(new_names.len(), self.names.len(), "new_names must have the same length as the current names");
        self.names = new_names;
    }

    pub fn serialize(&self, out: &mut impl Write) {
        self.tree.serialize(out).unwrap();
        for name in &self.names {
            let name_bytes = name.as_bytes();
            bincode::serialize_into(&mut *out, &(name_bytes.len() as u64)).unwrap();
            out.write_all(name_bytes).unwrap();
        }
    }

    pub fn load(input: &mut impl Read) -> Self {
        let tree = LcaTree::load(input).unwrap();
        let n = tree.n_nodes();
        let mut names = Vec::with_capacity(n);
        for _ in 0..n {
            let name_len: u64 = bincode::deserialize_from(&mut *input).unwrap();
            let mut name_bytes = vec![0_u8; name_len as usize];
            input.read_exact(&mut name_bytes).unwrap();
            names.push(String::from_utf8(name_bytes).unwrap());
        }
        Self { tree, names }
    }
}

const FEATURE_SET_FILE_MAGIC: [u8; 8] = *b"hksfs0.1";
const FEATURE_SET_FILE_VERSION: u32 = 1;

#[derive(Debug, Clone)]
pub struct Labeling<C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    pub color_assignments: C, // Map colex -> color id (integer)
    pub hierarchy: ColorHierarchy, // Color hierarchy for the color ids in color_assignments
    pub name: String, // Name of the feature set
}

impl<C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> Labeling<C> {
    fn serialize(&self, out: &mut impl Write) {
        self.color_assignments.serialize(out);
        self.hierarchy.serialize(out);

        bincode::serialize_into(&mut *out, &(self.name.as_bytes().len() as u64)).unwrap();
        out.write_all(self.name.as_bytes()).unwrap();
    }

    fn load(input: &mut impl Read) -> Self {
        let color_assignments = *C::load(input);
        let hierarchy = ColorHierarchy::load(input);

        let name_bytes_len: u64 = bincode::deserialize_from(&mut *input).unwrap();
        let mut name_bytes = vec![0_u8; name_bytes_len as usize];
        input.read_exact(&mut name_bytes).unwrap();
        let name = String::from_utf8(name_bytes).unwrap();

        Self { color_assignments, hierarchy, name }
    }

    pub fn serialize_to_file(&self, mut out: &mut impl Write) {
        out.write_all(&FEATURE_SET_FILE_MAGIC).unwrap();
        out.write_all(&FEATURE_SET_FILE_VERSION.to_le_bytes()).unwrap();
        self.serialize(out);
    }

    pub fn load_from_file(mut input: &mut impl Read) -> Self {
        let mut magic = [0_u8; 8];
        input.read_exact(&mut magic).unwrap();
        if magic != FEATURE_SET_FILE_MAGIC {
            panic!("Error loading feature set: invalid file format (magic constant mismatch)");
        }
        let mut version_bytes = [0_u8; 4];
        input.read_exact(&mut version_bytes).unwrap();
        let version = u32::from_le_bytes(version_bytes);
        if version != FEATURE_SET_FILE_VERSION {
            panic!("Error loading feature set: wrong version (found {}, expected {})", version, FEATURE_SET_FILE_VERSION);
        }
        Self::load(input)
    }
}

pub struct ColorStats {
    pub colored: usize,
    pub uncolored: usize,
    pub color_run_min: usize,
    pub color_run_max: usize,
    pub color_run_mean: f64,
    /// Count of SBWT positions assigned to each color node (indexed by color ID).
    pub color_counts: Vec<usize>,
}

pub struct SingleColoredKmersShort<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    inner: HksIndex<L,C>, // For each feature set, k-mers sharing an s-mer have been made to have the same color: the LCA in the color hierarchy
    k: usize, // the query k-mer length
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess + Sync + Send, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> SingleColoredKmersShort<L,C> {

    pub fn new(mut index: HksIndex<L, C>, k: usize, n_threads: usize) -> Self {
        assert!(k <= index.sbwt.k());
        if k < index.sbwt.k() {
            let fs = &mut index.feature_set;
            log::info!("Preprocessing colors for {}-mer queries for feature set: {}", k, fs.name);
            fs.color_assignments.substitute_lca_for_s_mer_ranges(k, fs.hierarchy.tree(), &index.lcs, n_threads);
        }

        Self { inner: index, k }
    }

    pub fn query_k(&self) -> usize {
        self.k
    }

    pub fn into_inner(self) -> HksIndex<L, C> {
        self.inner
    }

    pub fn inner(&self) -> &HksIndex<L, C> {
        &self.inner
    }
}

pub struct KmerLookupIteratorShort<'a, 'b, L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    // This iterator should be initialized so that the first k-1 MS values are skipped
    matching_stats_iter: MatchingStatisticsIterator<'a, 'b, SbwtIndex::<SubsetMatrix>, L>,
    index: &'a SingleColoredKmersShort<L, C>,
    query_pattern_length: usize,
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> Iterator for KmerLookupIteratorShort<'_, '_, L, C> {
    type Item = Option<usize>; // Color id of k-mer

    fn next(&mut self) -> Option<Self::Item> {
        let (len, range) = self.matching_stats_iter.next()?;

        if len >= self.query_pattern_length {
            // s-mer is found in the sbwt
            assert!(range.len() > 0);

            // Because we have preprocessed to color array, all positions in the colex
            // range of an s-mer have the same color. So we can pick any of those. Let's
            // pick the first one.
            let color = self.index.inner.get_color(range.start);

            Some(color) // This is a Some(Some(color)). Means that the iterator produced something.
        } else {
            Some(None) // Iterator not finished but the k-mer is not found -> no color
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::color_storage::SimpleColorStorage;

    fn sequential_substitute(colors: &mut SimpleColorStorage, s: usize, hierarchy: &LcaTree, lcs: &LcsWrapper) {
        let n = colors.len();
        let mut run_start = 0usize;
        for colex in 1..=n {
            let run_continues = colex < n && lcs.get_lcs(colex) >= s;
            if !run_continues {
                if colex - run_start > 1 {
                    let mut merged: Option<usize> = None;
                    for pos in run_start..colex {
                        merged = hierarchy.lca_options(merged, colors.get_color(pos));
                    }
                    for pos in run_start..colex {
                        colors.set_color(pos, merged);
                    }
                }
                run_start = colex;
            }
        }
    }

    fn run_parallel_vs_sequential_substitute(seqs: &[Vec<u8>], k: usize, s: usize, n_threads: usize) {
        let slices: Vec<&[u8]> = seqs.iter().map(|v| v.as_slice()).collect();
        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::<sbwt::BitPackedKmerSortingMem>::new()
            .k(k)
            .build_lcs(true)
            .run(sbwt::SliceSeqStream::new(&slices));
        let lcs = lcs.unwrap();

        let color_names = (0..seqs.len()).map(|i| format!("label{i}")).collect();
        let hierarchy = ColorHierarchy::new_star(color_names);

        let streams: Vec<sbwt::VecSeqStream> = seqs.iter()
            .map(|seq| sbwt::VecSeqStream::new(std::slice::from_ref(seq)))
            .collect();

        let index: HksIndex<LcsWrapper, SimpleColorStorage> =
            crate::build::build(sbwt, lcs, streams, 1, hierarchy, "feature_set_name", None);

        // Sequential reference: run the simple single-threaded loop
        let (_, lcs_seq, feature_set) = index.clone().into_parts();
        let (mut colors, hierarchy) = (feature_set.color_assignments.clone(), feature_set.hierarchy.clone());
        let lcs_wrapper_seq = LcsWrapper::from(lcs_seq);
        sequential_substitute(&mut colors, s, hierarchy.tree(), &lcs_wrapper_seq);

        // Parallel: run substitute_lca_for_s_mer_ranges on the same initial color storage
        let (sbwt2, lcs2, feature_set_2) = index.into_parts();
        let (mut colors2, hierarchy2) = (feature_set_2.color_assignments.clone(), feature_set_2.hierarchy.clone());
        let lcs_wrapper = LcsWrapper::from(lcs2);
        colors2.substitute_lca_for_s_mer_ranges(s, hierarchy2.tree(), &lcs_wrapper, n_threads);

        let n = sbwt2.n_sets();
        for i in 0..n {
            assert_eq!(
                colors.get_color(i), colors2.get_color(i),
                "mismatch at colex {i}"
            );
        }
    }

    #[test]
    fn test_parallel_substitute_lca_small() {
        run_parallel_vs_sequential_substitute(
            &[b"ACGTACGT".to_vec(), b"TGCATGCA".to_vec()],
            4, 3, 4,
        );
    }

    #[test]
    fn test_parallel_substitute_lca_large() {
        // Generate two 10k pseudorandom DNA sequences
        let bases = b"ACGT";
        let mut seq_a = Vec::with_capacity(10_000);
        let mut seq_b = Vec::with_capacity(10_000);
        let mut state = 12345_u64;
        for _ in 0..10_000 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            seq_a.push(bases[((state >> 33) & 3) as usize]);
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            seq_b.push(bases[((state >> 33) & 3) as usize]);
        }
        run_parallel_vs_sequential_substitute(&[seq_a, seq_b], 15, 10, 8);
    }
}

// Wrapper so that we can implement the foreing trait RandomAccessU32
#[derive(Debug, Clone)]
pub struct LcsWrapper {
    pub inner: LcsArray
}

impl ContractLeft for LcsWrapper {
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        self.inner.contract_left(I, target_len)
    }
}

impl From<LcsArray> for LcsWrapper {
    fn from(lcs: LcsArray) -> Self {
        Self {inner: lcs}
    }
}

impl LcsAccess for LcsWrapper {
    fn get_lcs(&self, colex: usize) -> usize {
        self.inner.access(colex)
    }
}

impl RandomAccessU32 for LcsWrapper {
    fn len(&self) -> usize {
        self.inner.len()
    }

    fn get(&self, idx: usize) -> u32 {
        self.inner.access(idx) as u32
    }
}

impl MySerialize for LcsWrapper {
    fn serialize(&self, out: &mut impl Write) {
        self.inner.serialize(out).unwrap();
    }

    fn load(input: &mut impl Read) -> Box<Self> {
        let inner = LcsArray::load(input).unwrap();
        Box::new(LcsWrapper { inner })
    }
}
