#![allow(non_snake_case, clippy::needless_range_loop, clippy::len_zero)] // Using upper-case variable names from the source material

use std::{collections::HashMap, fs::File, io::{BufRead, BufReader, BufWriter, Read, Write}, path::{Path, PathBuf}, sync::{Arc, Mutex}};
use clap::{Parser, Subcommand};
use io::{LazyFileSeqStream, SingleSeqStream};
use jseqio::{reader::DynamicFastXReader, record::Record};
use sbwt::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, LcsArray, SbwtIndex, SbwtIndexVariant, SubsetMatrix, write_sbwt_index_variant};
use single_colored_kmers::{ColorHierarchy, SingleColoredKmers};
use parallel_queries::OutputWriter;

use crate::{color_storage::SimpleColorStorage, parallel_queries::RunWriter, single_colored_kmers::{ColorStats, LcsWrapper, SingleColoredKmersShort}};

mod single_colored_kmers;
mod lca_tree;
mod lca_support;
mod io;
mod parallel_queries;
mod single_threaded_queries;
mod util;
mod wavelet_tree;
mod traits;
mod color_storage;

type FixedKColorIndex = SingleColoredKmers<LcsWrapper, SimpleColorStorage>;

enum ColorIndex { // For now just one variant, might add more later
    FixedK(FixedKColorIndex),
}

// If these names change, remember to also update the hardcoded mention in the
// help text of the --labels argument in the Build subcommand below.
// The duplication exists because Rust's concat!() only accepts literals, so
// we cannot build a compile-time string from this slice.
static RESERVED_COLOR_NAMES: &[&str] = &["none"];

const HKS_FILE_ID: [u8; 8] = *b"hks0.1.3";
const FIXED_INDEX_TYPE_ID: [u8; 4] = *b"fixd";
//const FLEXIBLE_INDEX_TYPE_ID: [u8; 4] = *b"flex";

impl ColorIndex {
    fn serialize(&self, out: &mut impl Write) {
        match self {
            ColorIndex::FixedK(index) => {
                out.write_all(&HKS_FILE_ID).unwrap();
                out.write_all(&FIXED_INDEX_TYPE_ID).unwrap();
                index.serialize(out);
            },
        }
    }

    fn load(input: &mut impl Read) -> Self {
        let mut file_id = [0_u8; 8];
        input.read_exact(&mut file_id).unwrap();
        assert_eq!(file_id, HKS_FILE_ID, "Invalid HKS file ID (outdated index file?)");

        let mut type_id = [0_u8; 4];
        input.read_exact(&mut type_id).unwrap();
        match type_id {
            FIXED_INDEX_TYPE_ID => {
                let index = ColorIndex::FixedK(FixedKColorIndex::load(input));
                log::info!("Loaded index with s = {}", index.k());
                index
            },
            _ => {
                panic!("Unknown index type ID in HKS file: {}", String::from_utf8_lossy(&type_id));
            }
        }
    }

    fn k(&self) -> usize {
        match self {
            ColorIndex::FixedK(index) => index.k(),
        }
    }

    fn is_flexible(&self) -> bool {
        match self {
            ColorIndex::FixedK(_) => false,
        }
    }

    fn color_names(&self) -> &[String] {
        match self {
            ColorIndex::FixedK(index) => index.color_names(),
        }
    }

    fn n_colors_in_hierarchy(&self) -> usize {
        match self {
            ColorIndex::FixedK(index) => index.n_colors_in_hierarchy(),
        }
    }

    fn root_id(&self) -> usize {
        match self {
            ColorIndex::FixedK(index) => index.color_hierarchy().root(),
        }
    }

    fn color_hierarchy(&self) -> &crate::lca_tree::LcaTree {
        match self {
            ColorIndex::FixedK(index) => index.color_hierarchy(),
        }
    }

    fn n_kmers(&self) -> usize {
        match self {
            ColorIndex::FixedK(index) => index.n_kmers(),
        }
    }

    fn rename_labels(&mut self, new_names: Vec<String>) {
        match self {
            ColorIndex::FixedK(index) => index.rename_labels(new_names),
        }
    }

    fn color_stats(&self) -> ColorStats {
        match self {
            ColorIndex::FixedK(index) => index.color_stats(),
        }
    }
}

// It's allowed for there to be names in the hierarchy that are not in the provided names.
// But every provided name must be in the hierarchy.
// Returns the tree and the names in id order.
fn read_hierarchy_file(path: &PathBuf, provided_names: &[String]) -> (crate::lca_tree::LcaTree, Vec<String>) {

    for name in provided_names.iter() {
        if RESERVED_COLOR_NAMES.contains(&name.as_str()) {
            panic!("Error: can not use \"{}\" as a label name because it is a reserved name", name);
        }
    }

    // Build map: label -> id
    let mut name_to_id = HashMap::<&str, usize>::new();
    for name in provided_names.iter() {
        name_to_id.insert(name, name_to_id.len());
    }

    let lines = read_all_lines(path);

    // Read edges as (child, parent) name pairs; one edge per line
    let mut edges = Vec::<(usize, usize)>::new();
    for (i, line) in lines.iter().enumerate() {
        if line.trim().is_empty() { continue; }
        let mut parts = line.split_whitespace();
        let child_name = parts.next().unwrap_or_else(|| panic!("Hierarchy file: missing child name on line {i}"));
        let parent_name = parts.next().unwrap_or_else(|| panic!("Hierarchy file: missing parent name on line {i}"));
        let child_id = name_to_id.get(child_name).copied().unwrap_or_else(|| {
            let new_id = name_to_id.len();
            name_to_id.insert(child_name, new_id);
            new_id
        });
        let parent_id = name_to_id.get(parent_name).copied().unwrap_or_else(|| {
            let new_id = name_to_id.len();
            name_to_id.insert(parent_name, new_id);
            new_id
        });
        edges.push((child_id, parent_id));
    }

    for name in provided_names {
        assert!(name_to_id.contains_key(name.as_str()), "Provided label {} not found in hierarchy", name);
    }

    let n_nodes = name_to_id.len();
    // Collect all names, including the new ones we might have found during parsing the tree.
    let mut all_names: Vec<String> = vec![String::new(); n_nodes];
    name_to_id.iter().for_each(|(name, id)| all_names[*id] = name.to_string());

    let tree = crate::lca_tree::LcaTree::new(n_nodes, edges)
        .unwrap_or_else(|e| panic!("Invalid hierarchy file {}: {e}", path.display()));

    (tree, all_names)
}

fn build_hierarchy(hierarchy_path: &Option<PathBuf>, provided_names: Vec<String>) -> ColorHierarchy {
    if let Some(path) = hierarchy_path {
        let (tree, all_names) = read_hierarchy_file(path, &provided_names);
        ColorHierarchy::with_tree(tree, all_names)
    } else {
        // This will check that "root" is not used as a label
        ColorHierarchy::new_star(provided_names)
    }
}

fn add_colors<T: sbwt::SeqStream + Send>(
    sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: LcsArray,
    individual_streams: Vec<T>,
    n_threads: usize,
    out_path: PathBuf,
    hierarchy: ColorHierarchy,
    nones_to_multiples: bool,
) {
    log::info!("Marking colors");
    let mut index = FixedKColorIndex::new(sbwt, lcs, individual_streams, n_threads, hierarchy);
    if nones_to_multiples {
        log::info!("Turning Nones into roots");
        index.turn_nones_to_roots();
    }

    let index = ColorIndex::FixedK(index);

    log::info!("Writing to {}", out_path.display());
    let mut out = BufWriter::new(File::create(out_path.clone())
        .unwrap_or_else(|e| panic!("Could not create output file {}: {e}", out_path.display())));
    index.serialize(&mut out);
    let out_size = std::fs::metadata(&out_path).unwrap().len() as f64;
    log::info!("Index size on disk: {}", human_bytes::human_bytes(out_size));
}

#[derive(Parser)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Subcommands,
}

#[derive(Subcommand)]
pub enum Subcommands {
    #[command(arg_required_else_help = true)]
    Build {
        #[arg(short, required = true, default_value = "31", help = "Maximum query length, up to 256. Warning: using a large value of s takes a lot of memory or disk during construction.", value_parser = clap::value_parser!(u64).range(1..=256))] // 256 is an upper limit of SBWT
        s: u64,

        #[arg(help = "A file with one fasta/fastq filename per line, one per label", long, help_heading = "Input", conflicts_with = "label_by_seq")]
        label_by_file: Option<PathBuf>,

        #[arg(help = "Give input as a single FASTA file, one sequence per label", long, help_heading = "Input", conflicts_with = "label_by_file")]
        label_by_seq: Option<PathBuf>,

        #[arg(help = "Optional: a fasta/fastq file containing the unitigs of all the k-mers in the input files. More generally, any sequence file with same k-mers will do (unitigs, matchtigs, eulertigs...). This speeds up construction and reduces the RAM and disk usage", short, long, help_heading = "Input")]
        unitigs: Option<PathBuf>,

        #[arg(help = "Output filename", short, long, required = true)]
        output: PathBuf,

        #[arg(help = "Run in external memory construction mode using the given directory as temporary working space. This reduces the RAM peak but is slower. The resulting index will still be exactly the same.", long = "external-memory")]
        temp_dir: Option<PathBuf>,

        #[arg(help = "Do not add reverse complemented k-mers", long = "forward-only")]
        forward_only: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4", value_parser = clap::value_parser!(u64).range(1..))]
        n_threads: u64,

        #[arg(help = "Optional: a precomputed Bit Matrix SBWT file of the input k-mers. Must have been built with --add-all-dummy-paths", long = "load-sbwt", help_heading = "Advanced use")]
        sbwt_path: Option<PathBuf>,

        #[arg(help = "Optional: a precomputed LCS file of the optional SBWT file. Must have been built with --add-all-dummy-paths", long = "load-lcs", help_heading = "Advanced use")]
        lcs_path: Option<PathBuf>,

        // The reserved names are hardcoded here because concat!() only accepts literals, not slice elements.
        // If RESERVED_COLOR_NAMES changes, update this help text accordingly.
        #[arg(help = "Optional: a file with one label name per line, in the same order as the input files. Defaults to using the input filenames as labels. The label \"none\" is reserved and cannot be used.", long = "labels", help_heading = "Input")]
        labels: Option<PathBuf>,

        #[arg(help = "Optional: a file describing the label hierarchy tree. Defaults to a star (all labels as children of a single root, named \"root\").", long = "hierarchy", help_heading = "Input")]
        hierarchy: Option<PathBuf>,

        #[arg(help = "Hidden option: After building, turn all \"none\" labels into \"multiple\"", long = "none-to-multiple", default_value = "false", hide = true)]
        none_to_multiple: bool,

        #[arg(help = "Optional: save the SBWT and LCS arrays to the given path prefix (writes <prefix>.sbwt and <prefix>.lcs).", long = "save-sbwt-and-lcs", help_heading = "Advanced use")]
        sbwt_and_lcs_save_prefix: Option<PathBuf>,

    },

    #[command(arg_required_else_help = true)]
    Lookup {
        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4", value_parser = clap::value_parser!(u64).range(1..))]
        n_threads: u64,

        #[arg(help = "Query k-mer length. Must be less or equal to the value of s used in index construction. If not given, defaults to the same k as during index construction.", short, required = false, value_parser = clap::value_parser!(u64).range(1..=256))] // 256 is an upper limit of SBWT
        k: Option<u64>,

        #[arg(help = "Print query names instead of query rank integers.", long = "report-query-names")]
        report_query_names: bool,

        #[arg(help = "Print lines for runs of k-mers not found in the index. The miss symbol is '-' normally, or 'none' when --report-label-names is set.", long = "report-misses")]
        report_misses: bool,

        #[arg(help = "Do not print the header line.", long = "no-header")]
        no_header: bool,

        #[arg(help = "Number of bases processed per batch in parallel query execution. Increasing this value increases RAM usage but may improve query time and/or parallelism.", long = "batch-size", default_value = "1000000", help_heading = "Advanced", value_parser = clap::value_parser!(u64).range(1..))]
        batch_size: u64,

        #[arg(help = "Report internal label id integers instead of label names. This might save a lot of space if the labels are long. Use --print-hierarchy to print the internal ids.", long = "report-label-ids", help_heading = "Advanced")]
        report_label_ids: bool,

    },

    #[command(about = "Print statistics about an index file.")]
    Stats {
        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,
    },

    #[command(about = "Print how the number of s-mers for each node in the hierarchy, for all 1 <= k <= s")]
    NodeStats {
        #[arg(help = "Path to the index file", long, required = true)]
        index: PathBuf,

        #[arg(help = "Print internal label ids instead of label names.", long = "report-label-ids")]
        report_color_ids: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,
    },

    #[command(about = "Print the label hierarchy of an index file. Output: number of labels on the first line, then all label names one per line (ids 0,1,2...), then all edges as space-separated child parent pairs, one per line.")]
    PrintHierarchy {
        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,
    },

    #[command(arg_required_else_help = true, about = "Simple reference implementation for debugging this program.")]
    LookupDebug {
        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,
    },

    #[command(arg_required_else_help = true, hide = true, about = "Rename labels in an index file.")]
    RenameLabels {
        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Path to a TSV file with two columns: label internal id (integer) and new label name", long = "new-names", required = true)]
        new_names: PathBuf,

        #[arg(help = "Output filename for the updated index. Can be the same as the input filename to update in place.", short, long, required = true)]
        output: PathBuf,
    },

}

struct DynamicFastXReaderWrapper {
    inner: DynamicFastXReader,
}

impl sbwt::SeqStream for DynamicFastXReaderWrapper{
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.inner.read_next().unwrap().map(|x| x.seq)
    }
}

fn load_seq_names(query_path: &PathBuf) -> Vec<String> {
    log::info!("Collecting sequence names from {} ...", query_path.display());
    let mut name_reader = DynamicFastXReader::from_file(query_path)
        .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));
    let mut seq_names = Vec::new();
    while let Some(rec) = name_reader.read_next().unwrap() {
        let header = std::str::from_utf8(rec.head).unwrap();
        let name = header.split_whitespace().next().unwrap_or(header);
        seq_names.push(name.to_string());
    }
    log::info!("Collected {} sequence names from query file", seq_names.len());
    seq_names
}

fn run_queries<W: RunWriter>(n_threads: usize, reader: DynamicFastXReader, index: ColorIndex, batch_size: usize, k: usize, writer: W) {
    let reader = DynamicFastXReaderWrapper { inner: reader };
    let ColorIndex::FixedK(index) = index;
    if k < index.k() {
        log::info!("Preprocessing colors for {}-mer queries", k);
        let s_index = SingleColoredKmersShort::new(index, k, n_threads);
        log::info!("Running {}-mer queries", k);
        parallel_queries::lookup_parallel(n_threads, reader, &s_index, batch_size, k, writer);
    } else {
        parallel_queries::lookup_parallel(n_threads, reader, &index, batch_size, k, writer);
    }
}

fn compute_node_stats(index: ColorIndex, report_color_names: bool, n_threads: usize) {
    use rayon::prelude::*;

    let color_names: Option<Vec<String>> = report_color_names.then(|| index.color_names().to_vec());
    let ColorIndex::FixedK(mut index) = index;
    let k = index.k();

    log::info!("Preprocessing: marking dummy nodes");
    let dummy_marks = index.sbwt().compute_dummy_node_marks();
    log::info!("Preprocessing: Building SBWT select support");
    index.build_sbwt_select();

    let stdout_mutex = std::sync::Mutex::new(std::io::BufWriter::new(std::io::stdout()));
    { let mut h = stdout_mutex.lock().unwrap(); writeln!(h, "k\tlabel\tcount").unwrap(); h.flush().unwrap(); }

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(|| {
        let k_values: Vec<usize> = (1..=k).rev().collect(); // Need to collect because par_iter does not take rev()
        k_values.into_par_iter().for_each(|s| {
            log::info!("Computing node stats for s = {}", s);
            let counts = index.node_stats(s, &dummy_marks);
            let mut out = String::new();
            for color in 0..counts.len() {
                let color_label = if let Some(ref names) = color_names {
                    names[color].clone()
                } else {
                    color.to_string()
                };
                out.push_str(&format!("{}\t{}\t{}\n", s, color_label, counts[color]));
            }
            let mut stdout = stdout_mutex.lock().unwrap();
            stdout.write_all(out.as_bytes()).unwrap();
            stdout.flush().unwrap();
        });
    });
}

fn rename_labels(index_path: &PathBuf, new_names_path: &PathBuf, out_path: &PathBuf) {
    let mut index_input = BufReader::new(File::open(index_path)
        .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));
    let mut index = ColorIndex::load(&mut index_input);

    // Read current names and apply overrides from TSV
    let mut names: Vec<String> = index.color_names().to_vec();
    let tsv_reader = BufReader::new(File::open(new_names_path)
        .unwrap_or_else(|e| panic!("Could not open new names file {}: {e}", new_names_path.display())));
    for (line_num, line) in tsv_reader.lines().enumerate() {
        let line = line.unwrap();
        if line.trim().is_empty() { continue; }
        let mut cols = line.splitn(2, '\t');
        let id_str = cols.next().unwrap_or_else(|| panic!("Line {}: missing label id", line_num + 1));
        let new_name = cols.next().unwrap_or_else(|| panic!("Line {}: missing new name", line_num + 1)).trim_end_matches(['\n', '\r']).to_owned();
        let id: usize = id_str.trim().parse().unwrap_or_else(|_| panic!("Line {}: label id is not a valid integer: {}", line_num + 1, id_str));
        if id >= names.len() {
            panic!("Line {}: label id {} is out of range (index has {} labels)", line_num + 1, id, names.len());
        }
        if RESERVED_COLOR_NAMES.contains(&new_name.as_str()) {
            panic!("Line {}: \"{}\" is a reserved label name and cannot be used", line_num + 1, new_name);
        }
        names[id] = new_name;
    }

    index.rename_labels(names);

    log::info!("Writing updated index to {}", out_path.display());
    let mut out = BufWriter::new(File::create(out_path)
        .unwrap_or_else(|e| panic!("Could not create output file {}: {e}", out_path.display())));
    index.serialize(&mut out);
}

// Load SBWT and LCS, or build from scratch if not given
fn save_sbwt_and_lcs_if_requested(sbwt: &SbwtIndexVariant, lcs: &LcsArray, prefix: &Option<PathBuf>) {
    
    if let Some(prefix) = prefix {
        let sbwt_out_path = PathBuf::from({ let mut s = prefix.as_os_str().to_os_string(); s.push(".sbwt"); s });
        let lcs_out_path = PathBuf::from({ let mut s = prefix.as_os_str().to_os_string(); s.push(".lcs"); s });
        if let Some(parent) = sbwt_out_path.parent() {
            std::fs::create_dir_all(parent).unwrap();
        }
        log::info!("Saving SBWT to {}", sbwt_out_path.display());
        let mut sbwt_out = BufWriter::new(File::create(&sbwt_out_path)
            .unwrap_or_else(|e| panic!("Could not create SBWT file {}: {e}", sbwt_out_path.display())));
        write_sbwt_index_variant(sbwt, &mut sbwt_out).unwrap();
        log::info!("Saving LCS to {}", lcs_out_path.display());
        let mut lcs_out = BufWriter::new(File::create(&lcs_out_path)
            .unwrap_or_else(|e| panic!("Could not create LCS file {}: {e}", lcs_out_path.display())));
        lcs.serialize(&mut lcs_out).unwrap();
    }
}

fn get_sbwt_and_lcs(sbwt_path: Option<PathBuf>, lcs_path: Option<PathBuf>, temp_dir: Option<PathBuf>, all_input_seqs: io::ChainedInputStream, n_threads: usize, add_rev_comps: bool, s: usize) -> (SbwtIndex<SubsetMatrix>, LcsArray){
    if let Some(sbwt_path) = sbwt_path {
        let mut input = BufReader::new(File::open(&sbwt_path)
            .unwrap_or_else(|e| panic!("Could not open SBWT file {}: {e}", sbwt_path.display())));
        let sbwt::SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt::load_sbwt_index_variant(&mut input).unwrap();
        log::info!("Loaded SBWT with {} k-mers", sbwt.n_kmers());
        let lcs = if let Some(lcs_path) = lcs_path {
            LcsArray::load(&mut BufReader::new(File::open(&lcs_path)
                .unwrap_or_else(|e| panic!("Could not open LCS file {}: {e}", lcs_path.display())))).unwrap()
        } else {
            LcsArray::from_sbwt(&sbwt, n_threads)
        };
        if sbwt.k() != s {
            panic!("The s specified ({}) does not match the k of the provided SBWT ({})", s, sbwt.k());
        }
        (sbwt, lcs)
    } else {
        let (sbwt, lcs) = if let Some(td) = temp_dir {
            // Use disk-based construction
            sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(add_rev_comps)
                .k(s)
                .build_lcs(true)
                .n_threads(n_threads)
                .precalc_length(8)
                .add_all_dummy_paths(true) // This is required for multi-k support
                .algorithm(BitPackedKmerSortingDisk::new().dedup_batches(false).temp_dir(&td))
            .run(all_input_seqs)
        } else {
            // Use in-memory construction
            sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(add_rev_comps)
                .k(s)
                .build_lcs(true)
                .n_threads(n_threads)
                .precalc_length(8)
                .add_all_dummy_paths(true) // This is required for multi-k support
                .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
            .run(all_input_seqs)
        };
        let lcs = lcs.unwrap(); // Ok because of build_lcs(true)
        (sbwt, lcs)
    }
}

fn read_all_lines(filename: &Path) -> Vec<String> {
    let reader = BufReader::new(File::open(filename)
        .unwrap_or_else(|e| panic!("Could not open input file {}: {e}", filename.display()))
    );
    
    let mut lines: Vec<String> = vec![];
    for line in reader.lines() {
        lines.push(line.unwrap())
    }
    lines
}

fn get_coloring_input_for_file_mode(file_of_files_path: &Path, labels_path: Option<&PathBuf>, hierarchy_path: &Option<PathBuf>, add_rev_comps: bool) -> (ColorHierarchy, Vec<LazyFileSeqStream>){
    let input_paths: Vec<PathBuf> = read_all_lines(file_of_files_path).into_iter().map(|f| PathBuf::from(f)).collect();

    // Read labels from file, or use filenames as default
    let labels: Vec<String> = if let Some(ref names_path) = labels_path {
        let names = read_all_lines(names_path);
        if names.len() != input_paths.len() {
            panic!("Label names file has {} names but there are {} input files", names.len(), input_paths.len());
        }
        names
    } else {
        // Use file paths as default labels
        read_all_lines(file_of_files_path)
    };
    let hierarchy = build_hierarchy(&hierarchy_path, labels);
    let individual_streams: Vec<LazyFileSeqStream> = input_paths.iter()
        .map(|p| LazyFileSeqStream::new(p.clone(), add_rev_comps))
        .collect();

    (hierarchy, individual_streams)
}

fn get_coloring_input_for_sequence_mode(seq_file: &Path, labels_path: Option<&PathBuf>, hierarchy_path: &Option<PathBuf>, add_rev_comps: bool) -> (ColorHierarchy, Vec<SingleSeqStream>) {
    // Read labels from file, or use sequence names as default
    let labels: Vec<String> = if let Some(ref names_path) = labels_path {
        log::info!("Reading label names from {}", names_path.display());
        read_all_lines(names_path)
    } else {
        log::info!("Reading sequence names from {}", seq_file.display());
        let mut pre_reader = DynamicFastXReader::from_file(&seq_file)
            .unwrap_or_else(|e| panic!("Could not open sequence file {}: {e}", seq_file.display()));
        let mut names = Vec::<String>::new();
        while let Some(rec) = pre_reader.read_next().unwrap() {
            names.push(String::from_utf8(rec.name().to_vec()).unwrap());
        }
        names
    };

    let n_labels = labels.len();
    let hierarchy = build_hierarchy(&hierarchy_path, labels);

    let shared_reader = Arc::new(Mutex::new(
        DynamicFastXReader::from_file(&seq_file)
            .unwrap_or_else(|e| panic!("Could not open sequence file {}: {e}", seq_file.display()))
    ));
    let individual_streams: Vec<io::SingleSeqStream> = (0..n_labels)
        .map(|_| io::SingleSeqStream::new(Arc::clone(&shared_reader), add_rev_comps))
        .collect();

    (hierarchy, individual_streams)

}

fn main() {

    if std::env::var("RUST_LOG").is_err() {
        // This is now unsafe since Rust 2024. Apparently
        // it's a flaw in Unix itself and cannot be called safely.
        unsafe { std::env::set_var("RUST_LOG", "info") }
    }
    env_logger::init();

    log::info!("Running hks version {}", env!("CARGO_PKG_VERSION"));

    let args = Cli::parse();

    match args.command {
        Subcommands::Build { label_by_file, label_by_seq, unitigs: unitigs_path, output: out_path, temp_dir, s, n_threads, forward_only, sbwt_path, lcs_path, labels: label_names_file, hierarchy: hierarchy_path, none_to_multiple, sbwt_and_lcs_save_prefix} => {

            let (s, n_threads) = (s as usize, n_threads as usize);

            // Create output directory if does not exist
            std::fs::create_dir_all(out_path.parent().unwrap()).unwrap();

            let add_rev_comps = !forward_only;

            // Determine SBWT inputs (all sequences together for k-mer set building)
            let sbwt_input_stream = if let Some(unitigs_path) = unitigs_path {
                io::ChainedInputStream::new(vec![unitigs_path.clone()])
            } else {
                let sbwt_input_paths: Vec<PathBuf> = if let Some(ref lbf) = label_by_file {
                    read_all_lines(lbf).into_iter().map(|line| PathBuf::from(line)).collect()
                } else {
                    vec![label_by_seq.as_ref().unwrap().clone()]
                };
                io::ChainedInputStream::new(sbwt_input_paths)
            };
            if let Some(fof) = label_by_file {
                // TODO: most of this code is duplicated in the else-branch. Refactor to extract the shared logic.
                // We load the coloring input first so we fail early if there is something wrong with it
                let (hierarchy, individual_streams) = get_coloring_input_for_file_mode(&fof, label_names_file.as_ref(), &hierarchy_path, add_rev_comps);
                let (sbwt, lcs) = get_sbwt_and_lcs(sbwt_path, lcs_path, temp_dir, sbwt_input_stream, n_threads, add_rev_comps, s);
                let sbwt_variant = SbwtIndexVariant::SubsetMatrix(sbwt); // Need to save in this form so that it has the type id like in sbwt-rs-cli
                save_sbwt_and_lcs_if_requested(&sbwt_variant, &lcs, &sbwt_and_lcs_save_prefix);
                let SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt_variant; // Get back the inner sbwt
                add_colors(sbwt, lcs, individual_streams, n_threads, out_path, hierarchy, none_to_multiple);
            } else {
                // We load the coloring input first so we fail early if there is something wrong with it
                let (hierarchy, individual_streams) = get_coloring_input_for_sequence_mode(&label_by_seq.unwrap(), label_names_file.as_ref(), &hierarchy_path, add_rev_comps);
                let (sbwt, lcs) = get_sbwt_and_lcs(sbwt_path, lcs_path, temp_dir, sbwt_input_stream, n_threads, add_rev_comps, s);
                let sbwt_variant = SbwtIndexVariant::SubsetMatrix(sbwt); // Need to save in this form so that it has the type id like in sbwt-rs-cli
                save_sbwt_and_lcs_if_requested(&sbwt_variant, &lcs, &sbwt_and_lcs_save_prefix);
                let SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt_variant; // Get back the inner sbwt
                add_colors(sbwt, lcs, individual_streams, n_threads, out_path, hierarchy, none_to_multiple);
            }

        },

        Subcommands::Lookup{query: query_path, index: index_path, n_threads, k, report_label_ids, report_query_names, report_misses, no_header, batch_size} => {

            let (n_threads, batch_size) = (n_threads as usize, batch_size as usize);
            
            log::info!("Loading the index ...");
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));

            let index_loading_start = std::time::Instant::now();
            let index = ColorIndex::load(&mut index_input);
            log::info!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());

            let k = k.unwrap_or(index.k() as u64) as usize;
            if k > index.k() {
                panic!("Error: query k = {} larger than indexing s = {}", k, index.k());
            }

            let color_names = if report_label_ids { None } else { Some(index.color_names().to_vec())};
            let seq_names = if report_query_names { Some(load_seq_names(&query_path)) } else { None };

            let reader = DynamicFastXReader::from_file(&query_path)
                .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));

            let stdout = BufWriter::with_capacity(1 << 21, std::io::stdout());
            let writer = OutputWriter::new(stdout, seq_names, color_names, report_misses, !no_header);

            log::info!("Running queries from {} ...", query_path.display());
            run_queries(n_threads, reader, index, batch_size, k, writer);
        },

        Subcommands::Stats { index: index_path } => {
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));
            let index = ColorIndex::load(&mut index_input);

            let stats = index.color_stats();
            println!("Index type:            {}", if index.is_flexible() { "flexible-k" } else { "fixed-k" });
            println!("k:                     {}", index.k());
            println!("Number of labels in hierarchy:      {}", index.n_colors_in_hierarchy());
            println!("Number of k-mers:      {}", index.n_kmers());
            println!("Labeled SBWT positions: {}", stats.colored);
            println!("Unlabeled SBWT positions:  {}", stats.uncolored);
            println!("Label run min length:  {}", stats.color_run_min);
            println!("Label run max length:  {}", stats.color_run_max);
            println!("Label run mean length: {:.2}", stats.color_run_mean);
            println!();
            println!("{:<10}  {}", "Count", "Label name");
            println!("{:<10}  {}", stats.uncolored, "none");
            for (id, name) in index.color_names().iter().enumerate() {
                println!("{:<10}  {}", stats.color_counts[id], name);
            }
        },

        Subcommands::NodeStats { index: index_path, report_color_ids, n_threads } => {
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));
            let index = ColorIndex::load(&mut index_input);
            compute_node_stats(index, !report_color_ids, n_threads);
        },

        Subcommands::PrintHierarchy { index: index_path } => {
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));
            let index = ColorIndex::load(&mut index_input);
            let names = index.color_names();
            let tree = index.color_hierarchy();
            let n = tree.n_nodes();
            println!("{}", n);
            for name in names {
                println!("{}", name);
            }
            for node in 0..n {
                if node != tree.root() {
                    println!("{} {}", names[node], names[tree.parent(node)]);
                }
            }
        },

        Subcommands::RenameLabels { index: index_path, new_names: new_names_path, output: out_path } => {
            rename_labels(&index_path, &new_names_path, &out_path);
        },

        Subcommands::LookupDebug{query: query_path, index: index_path} => {
            log::info!("Loading the index ...");
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));

            let index_loading_start = std::time::Instant::now();
            let index = ColorIndex::load(&mut index_input);
            log::info!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());
            log::info!("Running query debug implementation for {} ...", query_path.display());

            match index {
                ColorIndex::FixedK(index) => {
                    single_threaded_queries::lookup_single_threaded(&query_path, &index, index.k());
                },
            }
        }
    } 
}
