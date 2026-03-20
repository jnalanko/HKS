Some of these scripts use rust-script. Install with `cargo install rust-script`.

# Requirement:

- `sbwt` rust CLI installed
- `hks` installed
- `kraken2` installed

# Downloading the data

Everything except the reference and query genomes should be in the repository already.

The reference is the T2T assembly CHM13v2.0 and the query is the diploid T2T assembly hg002v1.1 Currently they can be downloaded like this: 


```bash
mkdir -p CHM13
curl https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -o CHM13/chm13v2.0.fa.gz

mkdir -p query
curl https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz -o query/hg002v1.1.fasta.gz
zcat query/hg002v1.1.fasta.gz > query/hg002v1.1.fasta
```

After this, run the following to split CHM13 into one file per chromosome and to extract information from fasta headers.

```bash
mkdir -p CHM13/chromosomes
mkdir -p feature_tree

# Split chromosomes to individual files
zcat CHM13/chm13v2.0.fa.gz | awk '
/^>/ {
    if (out) close(out)
    out = "CHM13/chromosomes/" n++ ".fna"
}
{ print > out }
'

zcat query/hg002v1.1.fasta.gz | grep ">" | sed 's/^>//' > query/names.txt
zcat query/hg002v1.1.fasta.gz > query/hg002v1.1.fasta
zcat CHM13/chm13v2.0.fa.gz | awk '/^>/ { print substr($1, 2) }' > feature_tree/chromosome_names.txt
```

Optional: verify sha1 checksums of all the data:

```bash
sha1sum query/hg002v1.1.fasta.gz query/names.txt feature_tree/chromosome_names.txt feature_tree/chromosome_tree.txt feature_tree/chromosome_fof.txt CHM13/chm13v2.0.fa.gz CHM13/chromosomes/0.fna CHM13/chromosomes/1.fna CHM13/chromosomes/2.fna CHM13/chromosomes/3.fna CHM13/chromosomes/4.fna CHM13/chromosomes/5.fna CHM13/chromosomes/6.fna CHM13/chromosomes/7.fna CHM13/chromosomes/8.fna CHM13/chromosomes/9.fna CHM13/chromosomes/10.fna CHM13/chromosomes/11.fna CHM13/chromosomes/12.fna CHM13/chromosomes/13.fna CHM13/chromosomes/14.fna CHM13/chromosomes/15.fna CHM13/chromosomes/16.fna CHM13/chromosomes/17.fna CHM13/chromosomes/18.fna CHM13/chromosomes/19.fna CHM13/chromosomes/20.fna CHM13/chromosomes/21.fna CHM13/chromosomes/22.fna CHM13/chromosomes/23.fna CHM13/chromosomes/24.fna > checksums/ours.txt

diff checksums/ours.txt checksums/correct.txt && echo "Checksums match"
```

# Pipeline 

```bash
mkdir -p temp
mkdir -p logs
mkdir -p index
mkdir -p out
mkdir -p plots
mkdir -p benchmark_results

# SBWT construction for HKS. Assumes we have 200 GB RAM available (-m 200)
mkdir CHM13/sbwt
/usr/bin/time -v sbwt build --input-list feature_tree/chromosome_fof.txt -t 64 -v -o CHM13/sbwt/CHM13-k63-added-dummies --temp-dir temp --in-memory --build-lcs -k 63 -m 200 -r --add-all-dummy-paths

python3 scripts/print_hks_build_benchmark_commands.py | bash
python3 scripts/print_hks_query_benchmark_commands.py | bash
python3 scripts/print_hks_to_heatmaps_commands.py | bash

# Kraken
python3 scripts/print_setup_kraken_dir_commands.py | bash
python3 scripts/print_kraken_minimizer_build_commands.py | bash
python3 scripts/print_kraken_minimizer_query_commands.py | bash
python3 scripts/print_kraken_to_heatmaps_commands.py | bash

# HKS parallel speedup
# Extract the first two chromosomes (chr1 maternal and paternal)
zcat query/hg002v1.1.fasta.gz | awk '/^>/{n++; if(n>2) exit} n' > query/hg002v1.1_chr1.fasta
# Run
python3 scripts/print_parallel_speedup_commands.py | bash

# Plots
hks print-hierarchy --index index/CHM13-s63.hks > feature_tree/hks_hierarchy_dump.txt
python3 scripts/draw_combined_heatmap.py
python3 scripts/plot_parallel_speedup.py
python3 scripts/draw_combined_heatmap.py preprint_results/heatmap_data/heatmap-k63-hks.tsv preprint_results/heatmap_data/heatmap-k63-m31-kraken.tsv plots/heatmap-combined.pdf

# Latex table
python3 scripts/build_latex_table.py
```

--

Detail: Kraken reserves id 0 to encode "not found" and id 1 for the root of the tree.

