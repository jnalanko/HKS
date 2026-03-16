# HKS: Hierarchical K-mer Sets

A data structure for **exact hierarchical variable-length k-mer annotation**. HKS assigns each k-mer in a collection to exactly one label in a user-defined category hierarchy, supports queries at any k-mer length k ≤ s from a single index, and provides comparable throughput to Kraken2 while being lossless.

Built on the [Spectral Burrows–Wheeler Transform (SBWT)](sbwt-rs-cli/api/README.md) for compact, high-speed DNA k-mer indexing.

> **Terminology note.** The manuscript uses *category* for the labels assigned to k-mers (e.g., chr1, LINE) and *category hierarchy* for the tree that organizes them. The CLI flags use the shorter synonym *color* (`--file-colors`, `--color-names`, `--report-color-names`) but the concepts are identical.

## Key Features

- **Variable-length queries**: Build once with `--max-k s`, then query at any k ≤ s — no need to rebuild for different k-mer lengths.
- **Exact, lossless annotation**: Every k-mer receives a precise label; no minimizer hashing, no Bloom filter approximations.
- **Hierarchical labeling**: A user-defined category hierarchy (e.g., chromosome morphology, repeat classification, taxonomy) resolves multi-matching k-mers to their lowest common ancestor (LCA).
- **Positional output**: Each position in the query receives a label, enabling detection of structural boundaries (e.g., translocations) within a single sequence.

## Installation

Requires a recent stable Rust toolchain. Clone with submodules, then build:

```bash
git clone --recurse-submodules https://github.com/jnalanko/HKS
cd HKS
cargo build --release
```

The binary is `target/release/hks`.

## Quick Start

```bash
# Build an index with max k-mer length 31, one category per input file
hks build \
  --max-k 31 \
  --file-colors example/file_of_files.txt \
  --output index.hks

# Query at the full k=31
hks lookup \
  --query example/query.fasta \
  --index index.hks

# Query the SAME index at k=15 — no rebuild needed
hks lookup \
  --query example/query.fasta \
  --index index.hks \
  -k 15
```

## Usage

### Build an index

Each input file (or each sequence, with `--sequence-colors`) defines one category. Provide a text file listing one FASTA/FASTQ path per line, one per category:

```bash
hks build \
  --max-k 63 \
  --file-colors input_files.txt \
  --color-names category_names.txt \
  --hierarchy hierarchy.txt \
  --output index.hks \
  --n-threads 32
```

| Flag | Description |
|------|-------------|
| `--max-k` (`-k`) | Maximum k-mer length (s). The index supports queries at any k ≤ s. |
| `--file-colors` | Text file with one FASTA/FASTQ path per category |
| `--sequence-colors` | FASTA/FASTQ where each sequence is one category |
| `--output` | Output index path |
| `--hierarchy` | Optional hierarchy file defining a tree over the categories |
| `--color-names` | Optional file with one category name per line |
| `--external-memory` | Directory for temporary files (reduces construction RAM usage) |
| `--forward-only` | Do not include reverse complements |
| `--n-threads` | Number of threads |

Input files can be gzip-compressed (`.gz`); decompression is handled automatically.

### Query k-mers

```bash
hks lookup \
  --query query.fasta \
  --index index.hks \
  -k 31 \
  --report-color-names \
  --report-query-names
```

| Flag | Description |
|------|-------------|
| `--query` | FASTA/FASTQ query file |
| `--index` | Index path |
| `-k` | Query k (must be ≤ max-k; defaults to max-k) |
| `--report-color-names` | Print category names instead of numeric IDs |
| `--report-query-names` | Include query sequence names in output |
| `--report-misses` | Also report positions with no k-mer match |
| `--no-header` | Suppress TSV header line |
| `--n-threads` | Number of threads |

Output is TSV with columns `query_rank`, `from_kmer`, `to_kmer`, `color`. Consecutive k-mers with the same label are run-length encoded into intervals.

### Smoothing

The raw output of `hks lookup` assigns each k-mer directly from the index. In practice, contiguous stretches of category-specific k-mers are often interrupted by short runs of multi-matching or novel k-mers (e.g., due to SNPs or indels between the query and the indexed sequences). The `hks smooth` subcommand applies a hierarchy-aware post-processing step that uses flanking context to recover more specific assignments.

The smoothing algorithm (Algorithm S1 in the paper) identifies windows exhibiting a specific → general → specific pattern in the hierarchy and reassigns interior intervals to the LCA of the flanking anchors. It iterates until convergence, typically within one or two passes.

```bash
# Pipe directly from lookup to smooth
hks lookup \
  --query assembly.fasta \
  --index index.hks \
  --report-color-names \
  --report-query-names \
  --report-misses \
| hks smooth --index index.hks \
> smoothed.tsv

# Or with files
hks smooth --index index.hks -i raw_output.tsv -o smoothed.tsv --max-gap 1000
```

| Flag | Description |
|------|-------------|
| `--index` | Path to the HKS index (used to load the category hierarchy) |
| `-i`, `--input` | Input TSV from `hks lookup` (default: stdin) |
| `-o`, `--output` | Output TSV file (default: stdout) |
| `--max-gap` | Maximum gap in bp between intervals for window detection (default: 1000) |

**Important:** The input to `hks smooth` must include unmatched positions. Always pass `--report-misses` to `hks lookup` when generating input for smoothing.

### Inspect an index

```bash
# Basic statistics
hks stats --index index.hks

# Per-node k-mer counts for all 1 ≤ s ≤ max-k (useful for choosing query k)
hks node-stats --index index.hks --report-color-names --n-threads 32
```

The `node-stats` output shows how many k-mers are assigned to each hierarchy node at each value of s. This helps determine where label specificity stabilizes (typically around k ≈ 20 for chromosome assignments).

### Hierarchy file format

By default, HKS uses a star topology: all categories are children of a single root node. The `--hierarchy` flag lets you supply a custom tree so that k-mers shared across multiple categories are reported at their Lowest Common Ancestor (LCA).

The file has three sections:

1. **Header line** – two integers: `<n_internal_nodes> <n_edges>`
2. **Internal node names** – one name per line, `n_internal_nodes` lines. These are the non-leaf nodes of the tree. IDs are assigned starting right after the leaf (category) IDs, in the order they appear here.
3. **Edge list** – one edge per line, `n_edges` lines. Each line is `<child_name> <parent_name>` (whitespace-separated). Both names must be either a category name (leaf) or an internal node name declared above.

Leaf names are the category names supplied via `--color-names`, or the input file paths when `--color-names` is omitted.

**Example** – three categories `A`, `B`, `C` grouped under two internal nodes:

```
Tree structure:

        root
       /    \
    clade1   C
    /    \
   A      B
```

Hierarchy file (assuming category names are `A`, `B`, `C`):

```
2 4
clade1
root
A clade1
B clade1
clade1 root
C root
```

## Reproducing Paper Results

The `chromosome_benchmark/` directory contains all scripts and configuration files needed to reproduce the results from the paper. See [`chromosome_benchmark/README.md`](chromosome_benchmark/README.md) for detailed instructions.

## Citing

If you use HKS in your work, please cite:

> Alanko JN, Ranallo-Benavidez TR, Barthel FP, Puglisi SJ, Marchet C. Hierarchical genomic feature annotation with variable-length queries. *RECOMB-Seq 2026*.

## License

[TODO: Add license information]
