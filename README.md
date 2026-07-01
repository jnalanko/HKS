# HKS: Hierarchical K-mer Sets

HKS is a variable-length k-mer index with hierarchical labels. The input consists of: 

* A label hierarchy described as a tree. For example, a phylogenetic tree.
* A set of k-mer sets, one for each label.

The index is built for a maximum s-mer length s, and allows queries for *any* k-mer length up to s. The query takes a sequence, and prints a file in bed-format annotating each input k-mer with the lowest common ancestor of the labels of that k-mer in the hierarchy.

Preprint available on [bioRXiv](https://www.biorxiv.org/content/10.64898/2026.03.15.711907v1). Note that this repository uses simplified terminology compared to the manuscript. Instead of the terms of "category" and "feature", this repository uses just "label" for all nodes in the label hierarchy. The next version of the manuscript will be updated to match the terminology in this repository.

## Installation

Requires a recent stable Rust toolchain. Clone with submodules, then build:

```bash
git clone --recurse-submodules https://github.com/jnalanko/HKS 
cd HKS
cargo build --release
```

The binary is `target/release/hks`.

## Usage

### Build an index

Index construction happens in two phases: first build the base k-mer index, then add a feature set labeling on top of it.

The `example/` directory contains a tiny example dataset with four files A.fna, B.fna, C.fna, D.fna. To index it, run:

```bash
# Phase 1: build the base index (all input k-mers, both strands)
hks build-base \
  -s 10 \
  --input-file-list example/file_of_files.txt \
  -o base.hksb

# Phase 2: add the feature set labeling
hks add-feature-set \
  -i base.hksb \
  -o features.hksf \
  --feature-file-list example/file_of_files.txt \
  --feature-hierarchy example/hierarchy.txt \
  --feature-set-name my_features \
  --variable-k-support
```

The `--variable-k-support` flag is required here because the query below uses a k-mer length (5) smaller than the index s (10). If you only ever query with k equal to s, you can omit this and get a more efficient index.

This creates two files: `base.hksb` (the base k-mer index) and `features.hksf` (the labeling of k-mers with the following hierarchy):

```
        root
       /    \
   clade2   D.fna
   /    \
 clade1  C.fna
 /    \
A.fna  B.fna
```

The full options for each phase are as follows:

```
Usage: hks build-base [OPTIONS] -s <S> --output <OUTPUT>

Options:
  -s <S>                            Maximum query length, up to 256. Warning: using a large value of s takes a lot of memory or disk during construction. [default: 31]
  -o, --output <OUTPUT>             Output filename. Recommended file extension: .hksb
      --external-memory <TEMP_DIR>  Run in external memory construction mode using the given directory as temporary working space. This reduces the RAM peak but is slower. The resulting index will still be exactly the same.
      --forward-only                Do not add reverse complemented k-mers
  -t, --n-threads <N_THREADS>       Number of parallel threads [default: 4]
      --mem-gigas <MEM_GIGAS>       RAM budget for SBWT construction in gigabytes. [default: 8]
  -h, --help                        Print help

Input:
      --input <INPUT>               Input fasta/fastq file. For multiple input files, see --input-file-list.
      --input-file-list <INPUT_FILE_LIST>
                                    A file with one input fasta/fastq filename per line.

Advanced use:
      --load-sbwt <SBWT_PATH>  Optional: a precomputed Bit Matrix SBWT file of the input k-mers. Must have been built with --add-all-dummy-paths
      --load-lcs <LCS_PATH>    Optional: a precomputed LCS file of the optional SBWT file. Must have been built with --add-all-dummy-paths
```

```
Usage: hks add-feature-set [OPTIONS] --index <INDEX> --output <OUTPUT> --feature-set-name <LABELING_NAME>

Options:
  -i, --index <INDEX>          Path to the existing base index file
  -o, --output <OUTPUT>        Output filename for the new feature set file
      --forward-only           Do not add reverse complemented k-mers
      --variable-k-support     Enable support for all k-mer lengths with k <= s in queries. Can not be used if feature priorities are given (--feature-priorities). This option requires that the base index has a dummy node representative for each prefix of the start of a sequence, otherwise the construction will crash with an error. Only use this if you are sure you know what you are doing.
  -t, --n-threads <N_THREADS>  Number of parallel threads [default: 4]
  -h, --help                   Print help

Features:
      --feature-file-list <LABEL_BY_FILE>
          A file with one fasta/fastq filename per line, one per feature. All k-mers in these files must already be present in the index.
      --feature-per-seq-file <LABEL_BY_SEQ>
          Give input as a single FASTA file, one sequence per feature. All k-mers in this file must already be present in the index.
      --feature-names <LABELS>
          Optional: a file with one feature name per line, in the same order as the input files/sequences. Defaults to using the input filenames or sequence names as features. The feature name "none" is reserved.
      --feature-hierarchy <HIERARCHY>
          Optional: a file describing the feature hierarchy tree. Defaults to a star (all features as children of a single root, named "root").
      --feature-set-name <LABELING_NAME>
          Name for the new feature set.
      --feature-priorities <NODE_PRIORITIES>
          Optional: a file assigning an integer priority to every node in the feature hierarchy (one "<name> <priority>" pair per line, whitespace-separated). Lower value = higher priority. Enables priority-aware LCA during construction. Nodes absent from the file default to priority 0.
```

### Query k-mers

To query the index built above with k-mer length 5 and the input file `example/query.fasta`, run the following:

```bash
hks lookup \
    -q example/query.fasta \
    -i base.hksb \
    --feature-set-file features.hksf \
    -k 5 \
    --report-query-names \
    --report-misses
```

This will print the following:

```
query_name	from_kmer	to_kmer	label_name
Q1	0	1	clade1
Q1	1	3	example/A.fna
Q1	3	4	clade1
Q1	4	5	example/B.fna
Q1	5	7	none
Q1	7	11	example/C.fna
Q1	11	12	clade2
Q1	12	15	root
Q2	0	1	clade1
Q2	1	8	example/B.fna
Q2	8	12	none
Q2	12	14	example/A.fna
Q2	14	15	clade2
Q2	15	19	none
Q2	19	22	example/C.fna
```

This means that k-mers `[0,1)` map to clade1, kmers `[1,3)` to A.fasta, kmers `[3,4)` to clade1 again, and so on.


The full query options are as follows:

```
Usage: hks lookup [OPTIONS] --index <INDEX> --query <QUERY>

Options:
  -i, --index <INDEX>
          Path to the base index file
      --feature-set-file <LABELING_FILE>
          Path to the feature set file. Defaults to the base index path with extension .hksf.
  -k <K>
          Query k-mer length. Must be less or equal to the value of s used in index construction. If not given, defaults to the same k as during index construction.
  -t, --n-threads <N_THREADS>
          Number of parallel threads [default: 4]
  -q, --query <QUERY>
          A fasta/fastq query file
      --report-query-names
          Print query names instead of query rank integers.
      --report-misses
          Print lines for runs of k-mers not found in the index. The miss symbol is 'none' normally, or '-' when --report-label-ids is set.
      --no-header
          Do not print the header line.
  -o, --output <OUTPUT>
          Output file. Defaults to stdout.
  -h, --help
          Print help

Advanced:
      --batch-size <BATCH_SIZE>  Number of bases processed per batch in parallel query execution. Increasing this value increases RAM usage but may improve query time and/or parallelism. [default: 1000000]
      --report-label-ids         Report internal label id integers instead of label names. This might save a lot of space if the labels are long. Use --print-hierarchy to print the internal ids.
```

### Hierarchy file format

By default, HKS uses a star topology: all features are children of a single root node. The `--feature-hierarchy` flag lets you supply a custom tree. The file is an edge list: one edge per line, each line is `<child feature> <parent feature>` (whitespace-separated). Every feature provided to the build command with `--feature-names` must appear in at least one edge.

See `example/hierarchy.txt` for an example.

### Other subcommands

Run the binary without arguments for more subcommands and their documentation.

### Citation

```
@article{alanko2026hierarchical,
    title        = {Hierarchical genomic feature annotation with variable-length queries},
    author       = {Alanko, Jarno N. and Ranallo-Benavidez, T. Rhyker and Barthel, Floris P. and Puglisi, Simon  J. and Marchet, Camille},
    year         = {2026},
    month        = {March},
    day          = {18},
    doi          = {10.64898/2026.03.15.711907},
    url          = {https://www.biorxiv.org/content/10.64898/2026.03.15.711907v1},
    publisher    = {bioRxiv},
    note         = {Preprint}
}
```
