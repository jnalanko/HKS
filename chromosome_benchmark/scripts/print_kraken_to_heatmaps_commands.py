print("set -xue")

for k in [15,31,47,63]:
    t = 16
    print(f"rust-script scripts/kraken_to_bed out/hg002-k{k}-t{t}-kraken.txt out/hg002-k{k}-kraken.bed")
    print(f"rust-script scripts/compute_heatmap_data.rs out/hg002-k{k}-kraken.bed out/heatmap-k{k}-kraken.tsv")
    print(f"python3 scripts/draw_heatmap.py out/heatmap-k{k}-kraken.tsv plots/heatmap-k{k}-kraken.png --kraken")
