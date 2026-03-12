# HKS — Disjoint K-mer Sets

A colored k-mer index for querying which input dataset (color) each k-mer originates from. Built on top of the [SBWT](sbwt-rs-cli/api/README.md) (Spectral Burrows-Wheeler Transform) for compact, high-speed DNA k-mer indexing.

## Overview

HKS assigns a **color** to each k-mer in an index.

- **Hierarchical colors** via a Lowest Common Ancestor (LCA) tree, allowing taxonomic or other organizational relationships between datasets.
- **Variable-k queries**: build at k, query at any s ≤ k.

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

Provide a text file listing one input FASTA/FASTQ path per line, one file per color.
The `example/` directory contains three single-color FASTA files and a file-of-files pointing to them:

```bash
hks build \
  --k 10 \
  --file-colors example/file_of_files.txt \
  --output index.hks
```

Key options:

| Flag | Description |
|------|-------------|
| `--k` | K-mer length |
| `--file-colors` | Text file with one FASTA/FASTQ path per color |
| `--sequence-colors` | FASTA/FASTQ where each sequence is one color |
| `--output` | Output index path |
| `--hierarchy` | Optional hierarchy file defining a tree with colors on the leaves |
| `--color-names` | Optional file with one color name per line |
| `--external-memory` | Directory for temporary files (reduces contruction RAM usage) |
| `--forward-only` | Do not include reverse complements |
| `--n-threads` | Number of threads |

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

Output is TSV with columns `query_rank`, `from_kmer`, `to_kmer`, `color` (run-length encoded intervals).

### Other subcommands

Run the binary without arguments for more subcommands and their documentation.
