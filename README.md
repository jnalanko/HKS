# HKS: Hierarchical K-mer Sets

HKS is a variable-length k-mer index with hierarchical labels. The input consists of: 

* A label hierarchy described as a tree. For example, a phylogenetic tree.
* A set of k-mer sets, one for each label.

The index is built for a maximum s-mer length s, and allows queries for *any* k-mer length up to s. The query takes a sequence, and prints a file in bed-format annotating each input k-mer with the lowest common ancestor of the labels of that k-mer in the hierarchy.

## Installation

Requires a recent stable Rust toolchain. Clone with submodules, then build:

```bash
git clone --recurse-submodules https://github.com/jnalanko/HKS 
cd HKS
cargo build --release
```

The binary is `target/release/hks`.

## Usage

### Build an index TODO update "color" -> "label"

The input to indexing is the maximum k-mer length s, and a file listing one input FASTA/FASTQ path per line, one file per color. Both DNA strands are indexed. By default, the names of the colors are the file paths of the input FASTA/FASTQ files. 

The `example/` directory contains a tiny example dataset with four files A.fna, B.fna, C.fna, D.fna. To index it, run:

```bash
hks build \
  -s 10 \
  --label-by-file example/file_of_files.txt \
  --hierarchy example/hierarchy.txt \
  --output index.hks
```

This indexes the data with the following color hierarchy.

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
Usage: hks build [OPTIONS] -s <S> --output <OUTPUT>

Options:
  -s <S>                            Maximum query length, up to 256. Warning: using a large value of s takes a lot of memory or disk during construction. [default: 31]
  -o, --output <OUTPUT>             Output filename
      --external-memory <TEMP_DIR>  Run in external memory construction mode using the given directory as temporary working space. This reduces the RAM peak but is slower. The resulting index will still be exactly the same.
      --forward-only                Do not add reverse complemented k-mers
  -t, --n-threads <N_THREADS>       Number of parallel threads [default: 4]
  -h, --help                        Print help

Input:
  -f, --file-colors <FILE_COLORS>
          A file with one fasta/fastq filename per line, one per color
      --sequence-colors <SEQUENCE_COLORS>
          Give input as a single file, one sequence per color
  -u, --unitigs <UNITIGS>
          Optional: a fasta/fastq file containing the unitigs of all the k-mers in the input files. More generally, any sequence file with same k-mers will do (unitigs, matchtigs, eulertigs...). This speeds up construction and reduces the RAM and disk usage
      --color-names <COLOR_NAMES_FILE>
          Optional: a file with one color name per line, in the same order as the input files. Defaults to using the input filenames as color names. The names "none" and "root" are reserved and cannot be used.
      --hierarchy <HIERARCHY>
          Optional: a file describing the color hierarchy tree. Defaults to a star (all colors as children of a single root, named "root").

Advanced use:
  -b, --sbwt-path <SBWT_PATH>  Optional: a precomputed Bit Matrix SBWT file of the input k-mers. Must have been built with --add-all-dummy-paths
  -l, --lcs-path <LCS_PATH>    Optional: a precomputed LCS file of the optional SBWT file. Must have been built with --add-all-dummy-paths
```

### Query k-mers TODO update "color" -> "label"

To query the index built above with k-mer length 5 and the input file `example/query.fasta`, run the following:

```bash
hks lookup \
    -q example/query.fasta \
    -i index.hks \
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

### Hierarchy file format

By default, HKS uses a star topology: all colors are children of a single root node. The `--hierarchy` flag lets you supply a custom tree so that k-mers are reported at their Lowest Common Ancestor (LCA) when they appear in multiple colors.

The file is an edge list: one edge per line, each line is `<child label> <parent label>` (whitespace-separated). Every label provided to the build command with `--labels` must appear in at least one edge. There can also be labels that were not provided with `--labels`.

See `example/hierarchy.txt` for an example.

### Other subcommands

Run the binary without arguments for more subcommands and their documentation.
