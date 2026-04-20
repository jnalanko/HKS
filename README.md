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

The input to indexing is the maximum k-mer length s, and a file listing one input FASTA/FASTQ path per line, one file per label. Both DNA strands are indexed. By default, the labels are the file paths of the input FASTA/FASTQ files. 

The `example/` directory contains a tiny example dataset with four files A.fna, B.fna, C.fna, D.fna. To index it, run:

```bash
hks build \
  -s 10 \
  --feature-file-list example/file_of_files.txt \
  --feature-hierarchy example/hierarchy.txt \
  --output-prefix index
```

This will create the index in two parts: `index.hksb` and `index.hksl`. The former contains
 the k-mer index, and the latter a labeling of the k-mers with the following hierarchy:

```
        root
       /    \
   clade2   D.fna
   /    \
 clade1  C.fna
 /    \
A.fna  B.fna
```

The full build options are as follows:

```
Usage: hks build [OPTIONS] -s <S> --output-prefix <OUTPUT_PREFIX>

Options:
  -s <S>                               Maximum query length, up to 256. Warning: using a large value of s takes a lot of memory or disk during construction. [default: 31]
  -o, --output-prefix <OUTPUT_PREFIX>  Output path prefix. Writes <PREFIX>.hksb (base index) and <PREFIX>.hksl (labeling).
      --external-memory <TEMP_DIR>     Run in external memory construction mode using the given directory as temporary working space. This reduces the RAM peak but is slower. The resulting index will still be exactly the same.
      --forward-only                   Do not add reverse complemented k-mers
  -t, --n-threads <N_THREADS>          Number of parallel threads [default: 4]
  -h, --help                           Print help

Features:
      --feature-file-list <LABEL_BY_FILE>
          A file with one fasta/fastq filename per line, one per feature
      --feature-per-seq-file <LABEL_BY_SEQ>
          Give input as a single FASTA file, one sequence per feature
      --feature-names <NAMES>
          Optional: a file with one feature name per line, in the same order as the input files/sequences. Defaults to using the input filenames or sequence names. The name "none" is reserved and cannot be used.
      --feature-hierarchy <HIERARCHY>
          Optional: a file describing the feature hierarchy tree. Defaults to a star (all features as children of a single root, named "root").
      --feature-priorities <NODE_PRIORITIES>
          Optional: a file assigning an integer priority to every node in the feature hierarchy (one "<name> <priority>" pair per line, whitespace-separated). Lower value = higher priority. Enables priority-aware LCA during construction, which keeps k-mers specific to high-priority subtrees rather than merging them to their common ancestor. Priorities are used during construction only and are not stored in the index. Warning: this makes construction use O(n^2) memory in the worst case, where n is the number of features in the hierarchy.
      --feature-set-name <LABELING_NAME>
          Name for the feature set [default: unnamed]

Advanced use:
  -u, --unitigs <UNITIGS>
          Optional: a fasta/fastq file containing the unitigs of all the k-mers in the input files. More generally, any sequence file with same k-mers will do (unitigs, matchtigs, eulertigs...). This speeds up construction and reduces the RAM and disk usage
      --load-sbwt <SBWT_PATH>
          Optional: a precomputed Bit Matrix SBWT file of the input k-mers. Must have been built with --add-all-dummy-paths
      --load-lcs <LCS_PATH>
          Optional: a precomputed LCS file of the optional SBWT file. Must have been built with --add-all-dummy-paths
      --save-sbwt-and-lcs <SBWT_AND_LCS_SAVE_PREFIX>
          Optional: save the SBWT and LCS arrays to the given path prefix (writes <prefix>.sbwt and <prefix>.lcs).
```

### Query k-mers

To query the index built above with k-mer length 5 and the input file `example/query.fasta`, run the following:

```bash
hks lookup \
    -q example/query.fasta \
    -i index.hksb \
    --labeling-file index.hksl \
    -k 5 \
    --report-query-names \
    --report-misses
```

This will print the following:

```
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
Usage: hks lookup [OPTIONS] --query <QUERY> --index <INDEX>

Options:
  -q, --query <QUERY>                          A fasta/fastq query file
  -i, --index <INDEX>                          Path to the base index file
      --labeling-file <LABELING_FILE>    Path to the labeling file
  -t, --n-threads <N_THREADS>                  Number of parallel threads [default: 4]
  -k <K>                                       Query k-mer length. Must be less or equal to the value of s used in index construction. If not given, defaults to the same k as during index construction.
      --report-query-names                     Print query names instead of query rank integers.
      --report-misses                          Print lines for runs of k-mers not found in the index. The miss symbol is 'none' normally, or '-' when --report-label-ids is set.
      --no-header                              Do not print the header line.
  -h, --help                                   Print help

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
