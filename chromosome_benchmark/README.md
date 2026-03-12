Some of these scripts use rust-script. Install with cargo install rust-script.

Everything except the reference and query genomes should be in the repository already.
The reference is CHM13 and the query is hg002. The CHM13 directory should look like this:

```
CHM13/
CHM13/sbwt
CHM13/sbwt/CHM13-k63-added-dummies.sbwt
CHM13/sbwt/CHM13-k63-added-dummies.lcs
CHM13/chromosomes
CHM13/chromosomes/7.fna
CHM13/chromosomes/8.fna
CHM13/chromosomes/23.fna
CHM13/chromosomes/9.fna
CHM13/chromosomes/0.fna
CHM13/chromosomes/19.fna
CHM13/chromosomes/22.fna
CHM13/chromosomes/16.fna
CHM13/chromosomes/6.fna
CHM13/chromosomes/11.fna
CHM13/chromosomes/21.fna
CHM13/chromosomes/14.fna
CHM13/chromosomes/5.fna
CHM13/chromosomes/18.fna
CHM13/chromosomes/10.fna
CHM13/chromosomes/24.fna
CHM13/chromosomes/13.fna
CHM13/chromosomes/15.fna
CHM13/chromosomes/17.fna
CHM13/chromosomes/4.fna
CHM13/chromosomes/20.fna
CHM13/chromosomes/1.fna
CHM13/chromosomes/12.fna
CHM13/chromosomes/2.fna
CHM13/chromosomes/3.fna
```

(The names of the chromosomes are at `feature_tree/chromosome_names.txt`)

The query directory should look like this:

```
query/
query/hg002v1.1.fasta.gz
query/names.txt
```

query/names.txt is just the fasta header of the query in order.

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

