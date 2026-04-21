#![allow(non_snake_case, clippy::needless_range_loop, clippy::len_zero)] // Using upper-case variable names from the source material

use std::{collections::HashMap, fs::File, io::{BufRead, BufReader, BufWriter, Read, Write}, path::{Path, PathBuf}, sync::{Arc, Mutex}};
use clap::{Parser, Subcommand};
use io::{LazyFileSeqStream, SingleSeqStream};
use jseqio::{reader::DynamicFastXReader, record::Record};
use sbwt::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, LcsArray, SbwtIndex, SbwtIndexVariant, SubsetMatrix, write_sbwt_index_variant};
use single_colored_kmers::{ColorHierarchy, Labeling, HksIndex};
use parallel_queries::OutputWriter;

use crate::{color_storage::SimpleColorStorage, lca_tree::LcaTree, parallel_queries::RunWriter, single_colored_kmers::{ColorStats, HksBase, LcsWrapper, SingleColoredKmersShort}, traits::ColoredKmerLookupAlgorithm};

mod single_colored_kmers;
mod build;
mod lca_tree;
mod lca_support;
mod priority_lca;
mod io;
mod parallel_queries;
mod single_threaded_queries;
mod util;
mod wavelet_tree;
mod traits;
mod color_storage;

type FixedKColorIndex = HksIndex<LcsWrapper, SimpleColorStorage>;
type ShortKColorIndex = SingleColoredKmersShort<LcsWrapper, SimpleColorStorage>;

enum ColorIndex { // For now just one variant, might add more later
    FixedK(FixedKColorIndex),
}

// If these names change, remember to also update the hardcoded mention in the
// help text of the --names argument in the Build subcommand below.
// The duplication exists because Rust's concat!() only accepts literals, so
// we cannot build a compile-time string from this slice.
static RESERVED_COLOR_NAMES: &[&str] = &["none"];

const HKS_FILE_ID: [u8; 8] = *b"hks0.1.5";
const FIXED_INDEX_TYPE_ID: [u8; 4] = *b"fixd";
//const FLEXIBLE_INDEX_TYPE_ID: [u8; 4] = *b"flex";

impl ColorIndex {
    fn serialize_base(&self, out: &mut impl Write) {
        match self {
            ColorIndex::FixedK(index) => {
                out.write_all(&HKS_FILE_ID).unwrap();
                out.write_all(&FIXED_INDEX_TYPE_ID).unwrap();
                index.base().serialize(out);
            },
        }
    }

    fn load(base_input: &mut impl Read, fs_input: &mut impl Read) -> Self {
        let mut file_id = [0_u8; 8];
        base_input.read_exact(&mut file_id).unwrap();
        assert_eq!(file_id, HKS_FILE_ID, "Invalid HKS file ID (outdated index file?)");

        let mut type_id = [0_u8; 4];
        base_input.read_exact(&mut type_id).unwrap();
        match type_id {
            FIXED_INDEX_TYPE_ID => {
                let base = HksBase::<LcsWrapper>::load(base_input);
                let labeling = Labeling::<SimpleColorStorage>::load_from_file(fs_input);
                let index = FixedKColorIndex::from_parts(base, labeling);
                log::info!("Loaded index with s = {}", index.k());
                ColorIndex::FixedK(index)
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

    fn n_kmers(&self) -> usize {
        match self {
            ColorIndex::FixedK(index) => index.n_kmers(),
        }
    }

    fn color_stats(&self) -> ColorStats {
        match self {
            ColorIndex::FixedK(index) => index.color_stats(),
        }
    }

    fn color_names(&self) -> &[String] {
        match self {
            ColorIndex::FixedK(index) => index.labeling().hierarchy.names(),
        }
    }

    fn color_hierarchy(&self) -> &LcaTree {
        match self {
            ColorIndex::FixedK(index) => index.labeling().hierarchy.tree(),
        }
    }

    fn n_colors_in_hierarchy(&self) -> usize {
        match self {
            ColorIndex::FixedK(index) => index.labeling().hierarchy.n_nodes(),
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

/// Parse a node-priority file of the form
///
/// ```text
/// node_name  priority
/// ```
///
/// one entry per line, tokens separated by whitespace. Returns a vector of
/// priorities indexed by node id (same ordering as `node_names`). Nodes that
/// do not appear in the file default to priority 0 (an INFO line is logged for
/// each such node). Unknown names and duplicate entries are errors.
fn parse_node_priorities(path: &Path, node_names: &[String]) -> Result<Vec<usize>, String> {
    let file = File::open(path).map_err(|e| format!("Could not open priorities file {}: {e}", path.display()))?;
    let name_to_id: HashMap<&str, usize> = node_names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();

    let mut priorities: Vec<Option<usize>> = vec![None; node_names.len()];
    for (lineno, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading {}:{}: {e}", path.display(), lineno + 1))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut toks = trimmed.split_whitespace();
        let name = toks.next().ok_or_else(|| format!("{}:{}: missing node name", path.display(), lineno + 1))?;
        let pri_tok = toks.next().ok_or_else(|| format!("{}:{}: missing priority for node {name}", path.display(), lineno + 1))?;
        if toks.next().is_some() {
            return Err(format!("{}:{}: expected exactly two tokens", path.display(), lineno + 1));
        }
        let pri: usize = pri_tok.parse().map_err(|e| format!("{}:{}: invalid priority {pri_tok:?}: {e}", path.display(), lineno + 1))?;
        let id = *name_to_id.get(name).ok_or_else(|| format!("{}:{}: unknown node name {name:?}", path.display(), lineno + 1))?;
        if priorities[id].is_some() {
            return Err(format!("{}:{}: duplicate priority for node {name:?}", path.display(), lineno + 1));
        }
        priorities[id] = Some(pri);
    }

    for (i, p) in priorities.iter_mut().enumerate() {
        if p.is_none() {
            log::info!("Node {:?} not found in priority file {}; defaulting to priority 0", node_names[i], path.display());
            *p = Some(0);
        }
    }

    Ok(priorities.into_iter().map(|p| p.unwrap()).collect())
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

fn resolve_labeling_file(index_path: &PathBuf, labeling_file: Option<PathBuf>) -> PathBuf {
    labeling_file.unwrap_or_else(|| index_path.with_extension("hksf"))
}

fn open_index(index_path: &PathBuf, labeling_file: Option<PathBuf>) -> ColorIndex {
    let labeling_path = resolve_labeling_file(index_path, labeling_file);
    let mut base_input = BufReader::new(File::open(index_path)
        .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));
    let mut fs_input = BufReader::new(File::open(&labeling_path)
        .unwrap_or_else(|e| panic!("Could not open labeling file {}: {e}", labeling_path.display())));
    ColorIndex::load(&mut base_input, &mut fs_input)
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
    BuildBase {
        #[arg(short, required = true, default_value = "31", help = "Maximum query length, up to 256. Warning: using a large value of s takes a lot of memory or disk during construction.", value_parser = clap::value_parser!(u64).range(1..=256))] // 256 is an upper limit of SBWT
        s: u64,

        #[arg(help = "Input fasta/fastq file. For multiple input files, see --input-file-list.", long, help_heading = "Input", conflicts_with = "input_file_list")]
        input: Option<PathBuf>,

        #[arg(help = "A file with one input fasta/fastq filename per line.", long, help_heading = "Input", conflicts_with = "input")]
        input_file_list: Option<PathBuf>,

        #[arg(help = "Output filename. Recommended file extension: .hksb", short = 'o', long = "output", required = true)]
        output: PathBuf,

        #[arg(help = "Run in external memory construction mode using the given directory as temporary working space. This reduces the RAM peak but is slower. The resulting index will still be exactly the same.", long = "external-memory")]
        temp_dir: Option<PathBuf>,

        #[arg(help = "Do not add reverse complemented k-mers", long = "forward-only")]
        forward_only: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4", value_parser = clap::value_parser!(u64).range(1..))]
        n_threads: u64,

        #[arg(help = "RAM budget for SBWT construction in gigabytes.", long = "mem-gigas", default_value = "8", value_parser = clap::value_parser!(u64).range(1..))]
        mem_gigas: u64,

        #[arg(help = "Optional: a precomputed Bit Matrix SBWT file of the input k-mers. Must have been built with --add-all-dummy-paths", long = "load-sbwt", help_heading = "Advanced use")]
        sbwt_path: Option<PathBuf>,

        #[arg(help = "Optional: a precomputed LCS file of the optional SBWT file. Must have been built with --add-all-dummy-paths", long = "load-lcs", help_heading = "Advanced use")]
        lcs_path: Option<PathBuf>,
    },

    AddFeatureSet {
        #[arg(help = "Path to the existing base index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Output filename for the new labeling file", short, long, required = true)]
        output: PathBuf,

        #[arg(help = "A file with one fasta/fastq filename per line, one per feature. All k-mers in these files must already be present in the index.", long = "feature-file-list", help_heading = "Features", conflicts_with = "label_by_seq")]
        label_by_file: Option<PathBuf>,

        #[arg(help = "Give input as a single FASTA file, one sequence per feature. All k-mers in this file must already be present in the index.", long = "feature-per-seq-file", help_heading = "Features", conflicts_with = "label_by_file")]
        label_by_seq: Option<PathBuf>,

        #[arg(help = "Optional: a file with one feature name per line, in the same order as the input files/sequences. Defaults to using the input filenames or sequence names as features. The feature name \"none\" is reserved.", long = "feature-names", help_heading = "Features")]
        labels: Option<PathBuf>,

        #[arg(help = "Optional: a file describing the feature hierarchy tree. Defaults to a star (all features as children of a single root, named \"root\").", long = "feature-hierarchy", help_heading = "Features")]
        hierarchy: Option<PathBuf>,

        #[arg(help = "Name for the new feature set.", long = "feature-set-name", required = true, help_heading = "Features")]
        labeling_name: String,

        #[arg(help = "Optional: a file assigning an integer priority to every node in the feature hierarchy (one \"<name> <priority>\" pair per line, whitespace-separated). Lower value = higher priority. Enables priority-aware LCA during construction. Nodes absent from the file default to priority 0.", long = "feature-priorities", help_heading = "Features")]
        node_priorities: Option<PathBuf>,

        #[arg(help = "Do not add reverse complemented k-mers", long = "forward-only")]
        forward_only: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4", value_parser = clap::value_parser!(u64).range(1..))]
        n_threads: u64,
    },


    #[command(arg_required_else_help = true)]
    Lookup {
        #[arg(help = "Path to the base index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Path to the labeling file. Defaults to the base index path with extension .hksf.", long = "labeling-file")]
        labeling_file: Option<PathBuf>,

        #[arg(help = "Query k-mer length. Must be less or equal to the value of s used in index construction. If not given, defaults to the same k as during index construction.", short, required = false, value_parser = clap::value_parser!(u64).range(1..=256))] // 256 is an upper limit of SBWT
        k: Option<u64>,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4", value_parser = clap::value_parser!(u64).range(1..))]
        n_threads: u64,

        #[command(flatten)] // These options are also used in the interactive prompt
        query_args: LookupQueryArgs,
    },

    // Hidden prompt command (hidden because the user interface might still change a lot)
    #[command(arg_required_else_help = true, hide = true, about = "Load an index once and run multiple queries interactively.")]
    Prompt {
        #[arg(help = "Path to the base index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Path to the labeling file. Defaults to the base index path with extension .hksf.", long = "labeling-file")]
        labeling_file: Option<PathBuf>,

        #[arg(help = "Query k-mer length for this session. Must be less or equal to the value of s used in index construction. If not given, defaults to the same k as during index construction.", short, required = false, value_parser = clap::value_parser!(u64).range(1..=256))]
        k: Option<u64>,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4", value_parser = clap::value_parser!(u64).range(1..))]
        n_threads: u64,
    },

    #[command(about = "Print statistics about an index file.")]
    Stats {
        #[arg(help = "Path to the base index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Path to the labeling file. Defaults to the base index path with extension .hksf.", long = "labeling-file")]
        labeling_file: Option<PathBuf>,
    },

    #[command(about = "Print how the number of s-mers for each node in the hierarchy, for all 1 <= k <= s")]
    NodeStats {
        #[arg(help = "Path to the base index file", long, required = true)]
        index: PathBuf,

        #[arg(help = "Path to the labeling file. Defaults to the base index path with extension .hksf.", long = "labeling-file")]
        labeling_file: Option<PathBuf>,

        #[arg(help = "Print internal label ids instead of label names.", long = "report-label-ids")]
        report_color_ids: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,
    },

    #[command(about = "Print the label hierarchy of an index file. Output: number of labels on the first line, then all label names one per line (ids 0,1,2...), then all edges as space-separated child parent pairs, one per line.")]
    PrintHierarchy {
        #[arg(help = "Path to the base index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Path to the labeling file. Defaults to the base index path with extension .hksf.", long = "labeling-file")]
        labeling_file: Option<PathBuf>,
    },

    #[command(arg_required_else_help = true, about = "Simple reference implementation for debugging this program.")]
    LookupDebug {
        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "Path to the base index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Path to the labeling file. Defaults to the base index path with extension .hksf.", long = "labeling-file")]
        labeling_file: Option<PathBuf>,
    },

}

#[derive(Parser, Debug)]
pub struct LookupQueryArgs {
    #[arg(help = "A fasta/fastq query file", short, long, required = true)]
    query: PathBuf,

    #[arg(help = "Print query names instead of query rank integers.", long = "report-query-names")]
    report_query_names: bool,

    #[arg(help = "Print lines for runs of k-mers not found in the index. The miss symbol is 'none' normally, or '-' when --report-label-ids is set.", long = "report-misses")]
    report_misses: bool,

    #[arg(help = "Do not print the header line.", long = "no-header")]
    no_header: bool,

    #[arg(help = "Number of bases processed per batch in parallel query execution. Increasing this value increases RAM usage but may improve query time and/or parallelism.", long = "batch-size", default_value = "1000000", help_heading = "Advanced", value_parser = clap::value_parser!(u64).range(1..))]
    batch_size: u64,

    #[arg(help = "Report internal label id integers instead of label names. This might save a lot of space if the labels are long. Use --print-hierarchy to print the internal ids.", long = "report-label-ids", help_heading = "Advanced")]
    report_label_ids: bool,

    #[arg(help = "Output file. Defaults to stdout.", short, long)]
    output: Option<PathBuf>,
}

struct DynamicFastXReaderWrapper {
    inner: DynamicFastXReader,
}

impl sbwt::SeqStream for DynamicFastXReaderWrapper{
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.inner.read_next().unwrap().map(|x| x.seq)
    }
}

fn open_fastx(path: &PathBuf) -> Result<DynamicFastXReader, String> {
    let file = File::open(path).map_err(|e| format!("Could not open file {}: {e}", path.display()))?;
    DynamicFastXReader::new(BufReader::new(file)).map_err(|e| format!("Could not read file {}: {e}", path.display()))
}

fn load_seq_names(query_path: &PathBuf) -> Result<Vec<String>, String> {
    log::info!("Collecting sequence names from {} ...", query_path.display());
    let mut name_reader = open_fastx(query_path)?;
    let mut seq_names = Vec::new();
    while let Some(rec) = name_reader.read_next().map_err(|e| format!("Error reading query file {}: {e}", query_path.display()))? {
        let header = std::str::from_utf8(rec.head).map_err(|e| format!("Invalid UTF-8 in sequence header in {}: {e}", query_path.display()))?;
        let name = header.split_whitespace().next().unwrap_or(header);
        seq_names.push(name.to_string());
    }
    log::info!("Collected {} sequence names from query file", seq_names.len());
    Ok(seq_names)
}

struct LookupAlgorithmImpl<'a> {
    index: &'a ShortKColorIndex,
}

impl<'a> ColoredKmerLookupAlgorithm for LookupAlgorithmImpl<'a> {
    fn lookup_kmers(&self, query: &[u8], k: usize) -> impl Iterator<Item = Option<usize>> {
        self.index.inner().lookup_kmers(query, k)
    }
}

fn run_queries<A: ColoredKmerLookupAlgorithm + Send + Sync, W: RunWriter>(n_threads: usize, reader: DynamicFastXReader, index: &A, batch_size: usize, k: usize, writer: W) {
    let reader = DynamicFastXReaderWrapper { inner: reader };
    parallel_queries::lookup_parallel(n_threads, reader, index, batch_size, k, writer);
}

fn run_lookup_with_args(index: &ShortKColorIndex, n_threads: usize, args: &LookupQueryArgs) -> Result<(), String> {
    let k = index.query_k();
    let seq_names = if args.report_query_names { Some(load_seq_names(&args.query)?) } else { None };
    let color_names: Option<Vec<String>> = if args.report_label_ids {
        None
    } else {
        Some(index.inner().labeling().hierarchy.names().to_vec())
    };
    let reader = open_fastx(&args.query)?;

    // A dynamic writer is fine performance-wise because it's wrapped in a buffered writer.
    let out: Box<dyn Write + Send> = if let Some(ref path) = args.output {
        Box::new(File::create(path).map_err(|e| format!("Could not create output file {}: {e}", path.display()))?)
    } else {
        Box::new(std::io::stdout())
    };
    let writer = OutputWriter::new(BufWriter::with_capacity(1 << 21, out), seq_names, color_names, args.report_misses, !args.no_header);

    let algo = LookupAlgorithmImpl { index };

    log::info!("Running queries from {} ...", args.query.display());
    run_queries(n_threads, reader, &algo, args.batch_size as usize, k, writer);
    Ok(())
}

fn run_prompt_loop(index: &ShortKColorIndex, n_threads: usize) {
    let stdin = std::io::stdin();
    let mut line = String::new();
    eprintln!("To run a query, type `-q example/query.fasta -o out.tsv` and press enter.");
    eprintln!("Hit enter without arguments for more instructions. ");
    loop {
        eprint!("lookup> ");
        line.clear();
        if stdin.lock().read_line(&mut line).unwrap() == 0 { break; }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            if let Err(e) = LookupQueryArgs::try_parse_from(["prompt", "--help"]) { eprintln!("{e}"); }
            continue;
        }
        if matches!(trimmed, "quit" | "exit" | "q") { break; }
        let tokens = std::iter::once("prompt").chain(trimmed.split_whitespace());
        match LookupQueryArgs::try_parse_from(tokens) {
            Ok(args) => if let Err(e) = run_lookup_with_args(index, n_threads, &args) { eprintln!("Error: {e}"); },
            Err(e) => eprintln!("{e}"),
        }
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

struct SbwtBuildOptions {
    seqs: io::ChainedInputStream,
    s: usize,
    add_rev_comps: bool,
    temp_dir: Option<PathBuf>,
    mem_gigas: usize,
}

enum SbwtSource {
    ComputeFromSeqs(SbwtBuildOptions),
    LoadFromDisk(PathBuf),
}

fn get_sbwt_and_lcs(sbwt_source: SbwtSource, lcs_path: Option<PathBuf>, n_threads: usize) -> (SbwtIndex<SubsetMatrix>, LcsArray){
    match sbwt_source {
        SbwtSource::ComputeFromSeqs(opts) => {
            // TODO: if lcs_path is given, use it instead of computing from the sbwt
            let (sbwt, lcs) = if let Some(td) = opts.temp_dir {
                // Use disk-based construction
                sbwt::SbwtIndexBuilder::new()
                    .add_rev_comp(opts.add_rev_comps)
                    .k(opts.s)
                    .build_lcs(true)
                    .n_threads(n_threads)
                    .precalc_length(8)
                    .add_all_dummy_paths(true) // This is required for multi-k support
                    .algorithm(BitPackedKmerSortingDisk::new().dedup_batches(false).temp_dir(&td).mem_gb(opts.mem_gigas))
                .run(opts.seqs)
            } else {
                // Use in-memory construction
                sbwt::SbwtIndexBuilder::new()
                    .add_rev_comp(opts.add_rev_comps)
                    .k(opts.s)
                    .build_lcs(true)
                    .n_threads(n_threads)
                    .precalc_length(8)
                    .add_all_dummy_paths(true) // This is required for multi-k support
                    .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false).mem_gb(opts.mem_gigas))
                .run(opts.seqs)
            };
            let lcs = lcs.unwrap(); // Ok because of build_lcs(true)
            (sbwt, lcs)
        },
        SbwtSource::LoadFromDisk(sbwt_path) => {
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
            (sbwt, lcs)
        }
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
        Subcommands::BuildBase { s, input, input_file_list, output: out_path, temp_dir, forward_only, n_threads, mem_gigas, sbwt_path, lcs_path  } => {

            let (s, n_threads, mem_gigas) = (s as usize, n_threads as usize, mem_gigas as usize);

            // Create output directory if does not exist
            if let Some(parent) = out_path.parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent).unwrap();
                }
            }

            // Open output file early to fail early if there is a problem
            let mut output_writer = BufWriter::new(File::create(&out_path).unwrap());

            let (sbwt, lcs) = match sbwt_path {
                Some(sbwt_path) => {
                    get_sbwt_and_lcs(SbwtSource::LoadFromDisk(sbwt_path), lcs_path, n_threads)
                }
                None => {
                    let add_rev_comps = !forward_only;

                    // Determine SBWT inputs (all sequences together for k-mer set building)
                    let sbwt_input_paths: Vec<PathBuf> = if let Some(ref input_fof) = input_file_list {
                        read_all_lines(input_fof).into_iter().map(PathBuf::from).collect()
                    } else {
                        vec![input.as_ref().unwrap().clone()]
                    };
                    let sbwt_input_stream = io::ChainedInputStream::new(sbwt_input_paths);
                    let sbwt_build_opts = SbwtBuildOptions{
                        seqs: sbwt_input_stream,
                        s,
                        add_rev_comps,
                        temp_dir,
                        mem_gigas,
                    };

                    get_sbwt_and_lcs(SbwtSource::ComputeFromSeqs(sbwt_build_opts), lcs_path, n_threads)
                }
            };

            // Package into HksBase
            let lcs = LcsWrapper::from(lcs);
            let base = HksBase::new(sbwt, lcs);

            log::info!("Writing base to {}", out_path.display());
            base.serialize(&mut output_writer);
        },

        Subcommands::Lookup { index: index_path, labeling_file, k, n_threads, query_args } => {
            log::info!("Loading the index ...");
            let index_loading_start = std::time::Instant::now();
            let index = open_index(&index_path, labeling_file);
            log::info!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());

            let ColorIndex::FixedK(index_inner) = index;
            let k = k.unwrap_or(index_inner.k() as u64) as usize;
            if k > index_inner.k() {
                panic!("Error: query k = {} larger than indexing s = {}", k, index_inner.k());
            }

            let n_threads = n_threads as usize;
            // Constructor does extra preprocessing if k < index_inner.k()
            let index = ShortKColorIndex::new(index_inner, k, n_threads);
            run_lookup_with_args(&index, n_threads, &query_args).unwrap_or_else(|e| panic!("{e}"));
        },

        Subcommands::Prompt { index: index_path, labeling_file, k, n_threads } => {
            log::info!("Loading the index ...");
            let index_loading_start = std::time::Instant::now();
            let index = open_index(&index_path, labeling_file);
            log::info!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());

            let ColorIndex::FixedK(index_inner) = index;
            let session_k = k.unwrap_or(index_inner.k() as u64) as usize;
            if session_k > index_inner.k() {
                panic!("Error: query k = {} larger than indexing s = {}", session_k, index_inner.k());
            }

            let n_threads = n_threads as usize;
            // Constructor does extra preprocessing if session_k < index_inner.k()
            let index = ShortKColorIndex::new(index_inner, session_k, n_threads);
            run_prompt_loop(&index, n_threads);
        },

        Subcommands::Stats { index: index_path, labeling_file } => {
            let index = open_index(&index_path, labeling_file);
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

        Subcommands::NodeStats { index: index_path, labeling_file, report_color_ids, n_threads } => {
            let index = open_index(&index_path, labeling_file);
            compute_node_stats(index, !report_color_ids, n_threads);
        },

        Subcommands::PrintHierarchy { index: index_path, labeling_file } => {
            let index = open_index(&index_path, labeling_file);
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

        Subcommands::AddFeatureSet { index: index_path, output: labeling_out_path, label_by_file, label_by_seq, labels: label_names_file, hierarchy: hierarchy_path, labeling_name, node_priorities: node_priorities_path, forward_only, n_threads } => {
            if label_by_file.is_none() && label_by_seq.is_none() {
                panic!("Error: one of --feature-file-list or --feature-per-seq-file is required");
            }

            let n_threads = n_threads as usize;
            let add_rev_comps = !forward_only;

            log::info!("Loading the base index ...");
            let mut base_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));
            let base = HksBase::<LcsWrapper>::load(&mut base_input);

            let labeling: Labeling<SimpleColorStorage> = if let Some(fof) = label_by_file {
                let (hierarchy, individual_streams) = get_coloring_input_for_file_mode(&fof, label_names_file.as_ref(), &hierarchy_path, add_rev_comps);
                let priorities = node_priorities_path.as_ref().map(|p| parse_node_priorities(p, hierarchy.names()).unwrap_or_else(|e| panic!("{e}")));
                build::build_labeling(&base, individual_streams, n_threads, hierarchy, &labeling_name, priorities)
            } else {
                let (hierarchy, individual_streams) = get_coloring_input_for_sequence_mode(&label_by_seq.unwrap(), label_names_file.as_ref(), &hierarchy_path, add_rev_comps);
                let priorities = node_priorities_path.as_ref().map(|p| parse_node_priorities(p, hierarchy.names()).unwrap_or_else(|e| panic!("{e}")));
                build::build_labeling(&base, individual_streams, n_threads, hierarchy, &labeling_name, priorities)
            };

            if let Some(parent) = labeling_out_path.parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent).unwrap();
                }
            }
            log::info!("Writing labeling to {}", labeling_out_path.display());
            let mut out = BufWriter::new(File::create(&labeling_out_path)
                .unwrap_or_else(|e| panic!("Could not create output file {}: {e}", labeling_out_path.display())));
            labeling.serialize_to_file(&mut out);
        },

        Subcommands::LookupDebug{query: query_path, index: index_path, labeling_file} => {
            log::info!("Loading the index ...");
            let index_loading_start = std::time::Instant::now();
            let index = open_index(&index_path, labeling_file);
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
