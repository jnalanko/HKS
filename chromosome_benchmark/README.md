Some of these scripts use rust-script. Install with `cargo install rust-script`.

# Downloading the data

Everything except the reference and query genomes should be in the repository already.

The reference is the T2T assembly CHM13v2.0 and the query is the diploid T2T assembly hg002v1.1 Currently they can be downloaded like this: 


```bash
mkdir -p CHM13
curl https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -o CHM13/chm13v2.0.fa.gz

mkdir -p query
curl https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz -o query/hg002v1.1.fasta.gz
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

zcat query/hg002v1.1.fasta.gz | grep ">" > query/names.txt
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

# SBWT construction for DKS. Assumes we have 200 GB RAM available (-m 200)
mkdir CHM13/sbwt
/usr/bin/time -v sbwt build --input-list feature_tree/chromosome_fof.txt -t 64 -v -o CHM13/sbwt/CHM13-k63-added-dummies --temp-dir temp --in-memory --build-lcs -k 63 -m 200 -r --add-all-dummy-paths

python3 scripts/print_dks_build_benchmark_commands.py | bash
python3 scripts/print_dks_query_benchmark_commands.py | bash
python3 scripts/print_dks_to_heatmaps_commands.py | bash

# Kraken
bash scripts/setup_kraken_dir.bash
python3 scripts/print_kraken_build_benchmark_commands.py | bash
python3 scripts/print_kraken_query_benchmark_commands.py | bash
python3 scripts/print_kraken_to_heatmaps_commands.py | bash

# DKS parallel speedup
seqtools extract-reads -r 0 -r 1 query/hg002v1.1.fasta.gz | gzip > query/hg002v1.1_chr1.fasta.gz
python3 scripts/print_parallel_speedup_commands.py | bash

# Plots
python3 scripts/draw_combined_heatmap.py
python3 scripts/plot_parallel_speedup.py
```

--

Detail: Kraken reserves id 0 to encode "not found" and id 1 for the root of the tree.

