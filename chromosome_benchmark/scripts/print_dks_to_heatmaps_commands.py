print("set -xue")

for k in [15,31,47,63]:
    t = 1
    print(f"rust-script scripts/compute_heatmap_data.rs out/hg002-k{k}-t{t}-dks.bed out/heatmap-k{k}-dks.tsv")
    print(f"python3 scripts/draw_heatmap.py out/heatmap-k{k}-dks.tsv plots/heatmap-k{k}-dks.png")
