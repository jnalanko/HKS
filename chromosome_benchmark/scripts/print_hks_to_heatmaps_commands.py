print("set -xue")

for k in [15,31,47,63]:
    t = 16
    print(f"rust-script scripts/compute_heatmap_data.rs out/hg002-k{k}-t{t}-hks.bed out/heatmap-k{k}-hks.tsv")
    print(f"python3 scripts/draw_heatmap.py out/heatmap-k{k}-hks.tsv plots/heatmap-k{k}-hks.png")
