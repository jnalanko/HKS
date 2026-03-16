# HKS: Hierarchical K-mer Sets

HKS is a variable-length k-mer index with hierarchical color labeling. The input consists of: 

* A set of k-mer sets, one for each color.
* A color hierarchy described as a tree where the colors are the leaves. For example, a phylogenetic tree.

The index is build for a maximum s-mer length s, and allows queries for *any* k-mer length up to s. The query takes a sequence, and prints a file in bed-format annotating each input k-mer with the lowest common ancestor of the colors of that k-mer in the hierarchy.

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

The input to indexing is the maximum k-mer length s, and a file listing one input FASTA/FASTQ path per line, one file per color. Both DNA strands are indexed. By default, the names of the colors are the file paths of the input FASTA/FASTQ files. 

The `example/` directory a tiny example dataset with three files A.fna, B.fna. C.fna. To index it, run:

```bash
hks build \
  -s 10 \
  --file-colors example/file_of_files.txt \
  --hierarchy example/hierarchy.txt \
  --output index.hks
```

This indexes the data with the following color hierarchy.

```
    root
   /    \
 clade1   C.fna
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

### Query k-mers

```bash
hks lookup \
  --query example/query.fasta \
  --index index.hks
```

Key options:

| Flag | Description |
|------|-------------|
| `--query` | FASTA/FASTQ query file |
| `--index` | Index path |
| `--k` | Query k (must be ≤ build k; defaults to build k) |
| `--report-color-names` | Print color names instead of IDs |
| `--report-query-names` | Include query sequence names in output |
| `--report-misses` | Also report positions with no k-mer match |
| `--no-header` | Suppress TSV header line |
| `--n-threads` | Number of threads |

Output is TSV with columns `query_rank`, `from_kmer`, `to_kmer`, `color`.

### Hierarchy file format

By default, HKS uses a star topology: all colors are children of a single root node. The `--hierarchy` flag lets you supply a custom tree so that k-mers are reported at their Lowest Common Ancestor (LCA) when they appear in multiple colors.

The file has three sections:

1. **Header line** – two integers: `<n_internal_nodes> <n_edges>`
2. **Internal node names** – one name per line, `n_internal_nodes` lines. These are the non-leaf nodes of the tree. IDs are assigned to them starting right after the leaf (color) IDs, in the order they appear here.
3. **Edge list** – one edge per line, `n_edges` lines. Each line is `<child_name> <parent_name>` (whitespace-separated). Both names must be either a color name (leaf) or an internal node name declared above.

Leaf names are the color names supplied via `--color-names`, or the input file paths when `--color-names` is omitted.

**Example** – three colors `A`, `B`, `C` grouped under two internal nodes:

```
Tree structure:

        root
       /    \
    clade1   C
    /    \
   A      B
```

Hierarchy file (assuming color names are `A`, `B`, `C`):

```
2 4
clade1
root
A clade1
B clade1
clade1 root
C root
```

### Other subcommands

Run the binary without arguments for more subcommands and their documentation.
