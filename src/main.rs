#![allow(non_snake_case, clippy::needless_range_loop, clippy::len_zero)] // Using upper-case variable names from the source material

use std::{collections::HashMap, fs::File, io::{BufRead, BufReader, BufWriter, Read, Write}, path::PathBuf, sync::{Arc, Mutex}};
use clap::{Parser, Subcommand};
use io::LazyFileSeqStream;
use jseqio::{reader::DynamicFastXReader, record::Record};
use sbwt::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, LcsArray};
use single_colored_kmers::{ColorHierarchy, SingleColoredKmers};
use parallel_queries::OutputWriter;

use crate::{color_storage::SimpleColorStorage, parallel_queries::RunWriter, single_colored_kmers::{ColorStats, LcsWrapper, SingleColoredKmersShort}};

mod single_colored_kmers;
mod lca_tree;
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
// help text of the --color-names argument in the Build subcommand below.
// The duplication exists because Rust's concat!() only accepts literals, so
// we cannot build a compile-time string from this slice.
static RESERVED_COLOR_NAMES: &[&str] = &["none", "root"];

const HKS_FILE_ID: [u8; 8] = *b"hks0.1.2";
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
        assert_eq!(file_id, HKS_FILE_ID, "Invalid HKS file ID");

        let mut type_id = [0_u8; 4];
        input.read_exact(&mut type_id).unwrap();
        match type_id {
            FIXED_INDEX_TYPE_ID => {
                let index = ColorIndex::FixedK(FixedKColorIndex::load(input));
                log::info!("Loaded index with k = {}", index.k());
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
}

// Returns the LcaTree and the internal node names (in node-ID order, after all leaf IDs).
fn read_hierarchy_file(path: &PathBuf, leaf_names: &[String]) -> (crate::lca_tree::LcaTree, Vec<String>) {

    for name in leaf_names.iter() {
        if RESERVED_COLOR_NAMES.contains(&name.as_str()) {
            panic!("Error: can not use \"{}\" as a color name because it is a reserved name", name);
        }
    }

    let mut name_to_id = HashMap::<&str, usize>::new();
    for (id, name) in leaf_names.iter().enumerate() {
        name_to_id.insert(name, id);
    }

    let file = File::open(path)
        .unwrap_or_else(|e| panic!("Could not open hierarchy file {}: {e}", path.display()));
    let mut lines = BufReader::new(file).lines();

    // First line: n_internal_nodes n_edges
    let first_line = lines.next()
        .unwrap_or_else(|| panic!("Hierarchy file {} is empty", path.display()))
        .unwrap();
    let mut parts = first_line.split_whitespace();
    let n_internal: usize = parts.next().unwrap().parse().unwrap();
    let n_edges: usize = parts.next().unwrap().parse().unwrap();

    // Read internal node names and assign IDs starting after leaf IDs
    let leaf_count = leaf_names.len();
    let mut internal_names = Vec::with_capacity(n_internal);
    for i in 0..n_internal {
        let name = lines.next()
            .unwrap_or_else(|| panic!("Hierarchy file: expected internal node name on line {}", i + 2))
            .unwrap();
        internal_names.push(name);
    }
    // Insert with stable string references from internal_names
    for (i, name) in internal_names.iter().enumerate() {
        name_to_id.insert(name.as_str(), leaf_count + i);
    }

    // Read edges
    let mut edges = Vec::with_capacity(n_edges);
    for i in 0..n_edges {
        let line = lines.next()
            .unwrap_or_else(|| panic!("Hierarchy file: expected edge on line {}", leaf_count + i + 2))
            .unwrap();
        let mut parts = line.split_whitespace();
        let child_name = parts.next().unwrap_or_else(|| panic!("Hierarchy file: missing child name in edge {i}"));
        let parent_name = parts.next().unwrap_or_else(|| panic!("Hierarchy file: missing parent name in edge {i}"));
        let child_id = *name_to_id.get(child_name)
            .unwrap_or_else(|| panic!("Hierarchy file: unknown node name '{child_name}'"));
        let parent_id = *name_to_id.get(parent_name)
            .unwrap_or_else(|| panic!("Hierarchy file: unknown node name '{parent_name}'"));
        edges.push((child_id, parent_id));
    }

    let n = leaf_count + n_internal;
    let tree = crate::lca_tree::LcaTree::new(n, edges)
        .unwrap_or_else(|e| panic!("Invalid hierarchy file {}: {e}", path.display()));
    (tree, internal_names)
}

fn build_hierarchy(hierarchy_path: &Option<PathBuf>, leaf_names: Vec<String>) -> ColorHierarchy {
    if let Some(path) = hierarchy_path {
        let (tree, internal_names) = read_hierarchy_file(path, &leaf_names);
        let mut all_names = leaf_names;
        all_names.extend(internal_names);
        ColorHierarchy::with_tree(tree, all_names)
    } else {
        ColorHierarchy::new_star(leaf_names)
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

        #[arg(help = "A file with one fasta/fastq filename per line, one per color", short, long, help_heading = "Input", conflicts_with = "sequence_colors")]
        file_colors: Option<PathBuf>,

        #[arg(help = "Give input as a single file, one sequence per color", long, help_heading = "Input", conflicts_with = "file_colors")]
        sequence_colors: Option<PathBuf>,

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

        #[arg(help = "Optional: a precomputed Bit Matrix SBWT file of the input k-mers. Must have been built with --add-all-dummy-paths", short = 'b', long, help_heading = "Advanced use")]
        sbwt_path: Option<PathBuf>,

        #[arg(help = "Optional: a precomputed LCS file of the optional SBWT file. Must have been built with --add-all-dummy-paths", short, long, help_heading = "Advanced use")]
        lcs_path: Option<PathBuf>,

        // The reserved names are hardcoded here because concat!() only accepts literals, not slice elements.
        // If RESERVED_COLOR_NAMES changes, update this help text accordingly.
        #[arg(help = "Optional: a file with one color name per line, in the same order as the input files. Defaults to using the input filenames as color names. The names \"none\" and \"root\" are reserved and cannot be used.", long = "color-names", help_heading = "Input")]
        color_names_file: Option<PathBuf>,

        #[arg(help = "Optional: a file describing the color hierarchy tree. Defaults to a star (all colors as children of a single root, named \"root\").", long = "hierarchy", help_heading = "Input")]
        hierarchy: Option<PathBuf>,

        #[arg(help = "Hidden option: After building, turn all \"none\" colors into \"multiple\"", long = "none-to-multiple", default_value = "false", hide = true)]
        none_to_multiple: bool,

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

        #[arg(help = "Print color names instead of color rank integers. K-mers present in multiple colors 'root' when this flag is set.", long = "report-color-names")]
        report_color_names: bool,

        #[arg(help = "Print query names instead of query rank integers.", long = "report-query-names")]
        report_query_names: bool,

        #[arg(help = "Print lines for runs of k-mers not found in the index. The miss symbol is '-' normally, or 'none' when --report-color-names is set.", long = "report-misses")]
        report_misses: bool,

        #[arg(help = "Do not print the header line.", long = "no-header")]
        no_header: bool,

        #[arg(help = "Number of bases processed per batch in parallel query execution. Increasing this value increases RAM usage but may improve query time and/or parallelism.", long = "batch-size", default_value = "1000000", help_heading = "Advanced", value_parser = clap::value_parser!(u64).range(1..))]
        batch_size: u64,
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

        #[arg(help = "Print color names instead of color ids", long = "report-color-names")]
        report_color_names: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,
    },

    #[command(arg_required_else_help = true, about = "Simple reference implementation for debugging this program.")]
    LookupDebug {
        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,
    },

    /* Outdated commmand. This was before the changes that introduced the color hierarchy
    #[command(arg_required_else_help = true, about = "Debug: build individual SBWTs per color and query them separately.")]
    IndividualSbwtDebug {
        #[arg(help = "A file with one fasta/fastq filename per line (one per color)", short, long, required = true, help_heading = "Input")]
        input: PathBuf,

        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "K-mer length (used for both indexing and querying)", short, required = true)]
        k: usize,

        #[arg(help = "Do not add reverse complemented k-mers", short = 'f', long = "forward-only")]
        forward_only: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,
    },
    */
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
    { let mut h = stdout_mutex.lock().unwrap(); writeln!(h, "s\tcolor\tcount").unwrap(); h.flush().unwrap(); }

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

// Reads a color names file with one name per line.
fn read_color_names_file(path: &PathBuf) -> Vec<String> {
    BufReader::new(File::open(path)
        .unwrap_or_else(|e| panic!("Could not open color names file {}: {e}", path.display())))
        .lines()
        .map(|l| l.unwrap())
        .collect()
}

/* This code is outdated. It is before the changes that introduced the color hierarchy tree 
fn individual_sbwt_debug(input_fof: &PathBuf, query_path: &PathBuf, k: usize, forward_only: bool, n_threads: usize) {
    // Read input file-of-files
    let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof)
        .unwrap_or_else(|e| panic!("Could not open input file {}: {e}", input_fof.display())))
        .lines().map(|f| PathBuf::from(f.unwrap())).collect();

    // Build one SBWT per input file, keeping all in memory
    let mut indices = Vec::new();
    for input_path in &input_paths {
        log::info!("Building SBWT for {}", input_path.display());
        let stream = LazyFileSeqStream::new(input_path.clone(), false);
        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(!forward_only)
            .k(k)
            .build_lcs(true)
            .n_threads(n_threads)
            .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
            .run(stream);
        indices.push((sbwt, lcs.unwrap()));
    }

    // Open query file
    let mut reader = DynamicFastXReader::from_file(query_path)
        .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));

    // Create OutputWriter (same format as lookup command)
    let stdout = BufWriter::with_capacity(1 << 17, std::io::stdout());
    let mut writer = OutputWriter::new(stdout, None, None, false, true);
    writer.write_header();

    // Stream query sequences, querying all SBWTs in lockstep
    let mut seq_id: isize = 0;
    while let Some(rec) = reader.read_next().unwrap() {
        let mut ms_iters: Vec<_> = indices.iter()
            .map(|(sbwt, lcs)| {
                let si = sbwt::StreamingIndex::new(sbwt, lcs);
                si.matching_statistics_iter(rec.seq)
            })
            .collect();

        // Skip first k-1 positions (not full k-mers yet)
        for _ in 0..k.saturating_sub(1) {
            for iter in ms_iters.iter_mut() { iter.next(); }
        }

        // Advance all iterators in lockstep, run-length encoding on the fly
        let mut run_start = 0usize;
        let mut run_color = None;
        let mut kmer_count = 0usize;
        loop {
            let steps: Vec<Option<_>> = ms_iters.iter_mut().map(|it| it.next()).collect();
            if steps.iter().any(|s| s.is_none()) { break; }

            // Union colors: a k-mer is in color i iff its MS length == k
            let color = steps.into_iter().enumerate().fold(
                None,
                |acc, (color_id, step)| {
                    let (len, _) = step.unwrap();
                    if len == k { acc.union(ColorVecValue::Single(color_id)) }
                    else { acc }
                },
            );

            if kmer_count == 0 {
                run_start = 0;
                run_color = color;
            } else if color != run_color {
                writer.write_run(seq_id, run_color, run_start..kmer_count);
                run_start = kmer_count;
                run_color = color;
            }
            kmer_count += 1;
        }
        if kmer_count > 0 {
            writer.write_run(seq_id, run_color, run_start..kmer_count);
        }
        seq_id += 1;
    }
    writer.flush();
}
*/

fn main() {

    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    log::info!("Running hks version {}", env!("CARGO_PKG_VERSION"));

    let args = Cli::parse();

    match args.command {
        Subcommands::Build { file_colors, sequence_colors, unitigs: unitigs_path, output: out_path, temp_dir, s, n_threads, forward_only, sbwt_path, lcs_path, color_names_file, hierarchy: hierarchy_path, none_to_multiple} => {

            let (s, n_threads) = (s as usize, n_threads as usize);

            // Create output directory if does not exist
            std::fs::create_dir_all(out_path.parent().unwrap()).unwrap();

            let add_rev_comps = !forward_only;

            // Determine SBWT input paths (all sequences together for k-mer set building)
            let sbwt_input_paths: Vec<PathBuf> = if let Some(ref fc) = file_colors {
                BufReader::new(File::open(fc)
                    .unwrap_or_else(|e| panic!("Could not open input file {}: {e}", fc.display())))
                    .lines().map(|l| PathBuf::from(l.unwrap())).collect()
            } else {
                vec![sequence_colors.as_ref().unwrap().clone()]
            };

            let (sbwt, lcs) = if let Some(sbwt_path) = sbwt_path {
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
                let all_input_seqs = if let Some(unitigs_path) = unitigs_path {
                    io::ChainedInputStream::new(vec![unitigs_path.clone()])
                } else {
                    io::ChainedInputStream::new(sbwt_input_paths)
                };
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
            };

            if let Some(fc) = file_colors {
                let input_paths: Vec<PathBuf> = BufReader::new(File::open(&fc)
                    .unwrap_or_else(|e| panic!("Could not open input file {}: {e}", fc.display())))
                    .lines().map(|l| PathBuf::from(l.unwrap())).collect();
                let color_names: Vec<String> = if let Some(ref names_path) = color_names_file {
                    let names = read_color_names_file(names_path);
                    if names.len() != input_paths.len() {
                        panic!("Color names file has {} names but there are {} input files", names.len(), input_paths.len());
                    }
                    names
                } else {
                    input_paths.iter().map(|p| p.as_os_str().to_str().unwrap().to_owned()).collect()
                };
                let hierarchy = build_hierarchy(&hierarchy_path, color_names);
                let individual_streams: Vec<LazyFileSeqStream> = input_paths.iter()
                    .map(|p| LazyFileSeqStream::new(p.clone(), add_rev_comps))
                    .collect();
                add_colors(sbwt, lcs, individual_streams, n_threads, out_path, hierarchy, none_to_multiple);
            } else {
                let sc = sequence_colors.unwrap();

                let color_names: Vec<String> = if let Some(ref names_path) = color_names_file {
                    log::info!("Reading color names from {}", names_path.display());
                    read_color_names_file(names_path)
                } else {
                    log::info!("Reading sequence names from {}", sc.display());
                    let mut pre_reader = DynamicFastXReader::from_file(&sc)
                        .unwrap_or_else(|e| panic!("Could not open sequence-colors file {}: {e}", sc.display()));
                    let mut names = Vec::<String>::new();
                    while let Some(rec) = pre_reader.read_next().unwrap() {
                        names.push(String::from_utf8(rec.name().to_vec()).unwrap());
                    }
                    names
                };

                let n_colors = color_names.len();
                let hierarchy = build_hierarchy(&hierarchy_path, color_names);

                let shared_reader = Arc::new(Mutex::new(
                    DynamicFastXReader::from_file(&sc)
                        .unwrap_or_else(|e| panic!("Could not open sequence-colors file {}: {e}", sc.display()))
                ));
                let individual_streams: Vec<io::SingleSeqStream> = (0..n_colors)
                    .map(|_| io::SingleSeqStream::new(Arc::clone(&shared_reader), add_rev_comps))
                    .collect();

                add_colors(sbwt, lcs, individual_streams, n_threads, out_path, hierarchy, none_to_multiple);
            }

        },

        Subcommands::Lookup{query: query_path, index: index_path, n_threads, k, report_color_names, report_query_names, report_misses, no_header, batch_size} => {

            let (n_threads, batch_size) = (n_threads as usize, batch_size as usize);
            
            log::info!("Loading the index ...");
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));

            let index_loading_start = std::time::Instant::now();
            let index = ColorIndex::load(&mut index_input);
            log::info!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());

            let k = k.unwrap_or(index.k() as u64) as usize;
            if k > index.k() {
                panic!("Error: query k = {} larger than indexing k = {}", k, index.k());
            }

            let color_names = report_color_names.then(|| index.color_names().to_vec());
            let seq_names = report_query_names.then(|| load_seq_names(&query_path));

            let reader = DynamicFastXReader::from_file(&query_path)
                .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));

            let stdout = BufWriter::with_capacity(1 << 21, std::io::stdout());
            let writer = OutputWriter::new(stdout, seq_names, color_names, index.root_id(), report_misses, !no_header);

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
            println!("Number of colors in hierarchy:      {}", index.n_colors_in_hierarchy());
            println!("Number of k-mers:      {}", index.n_kmers());
            println!("Colored SBWT positions: {}", stats.colored);
            println!("Uncolored SBWT positions:  {}", stats.uncolored);
            println!("Color run min length:  {}", stats.color_run_min);
            println!("Color run max length:  {}", stats.color_run_max);
            println!("Color run mean length: {:.2}", stats.color_run_mean);
            println!();
            println!("{:<10}  {}", "Count", "Color name");
            println!("{:<10}  {}", stats.uncolored, "none");
            for (id, name) in index.color_names().iter().enumerate() {
                println!("{:<10}  {}", stats.color_counts[id], name);
            }
        },

        Subcommands::NodeStats { index: index_path, report_color_names, n_threads } => {
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));
            let index = ColorIndex::load(&mut index_input);
            compute_node_stats(index, report_color_names, n_threads);
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
