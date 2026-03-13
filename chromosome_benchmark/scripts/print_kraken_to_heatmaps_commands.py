print("set -xue")

for k in [15,31,47,63]:
    for m in [15,22,32]:
        if m > k: continue
        t = 32
        print(f"rust-script scripts/kraken_to_bed out/hg002-k{k}-m{m}-t{t}-kraken.txt out/hg002-k{k}-m{m}-kraken.bed")
        print(f"rust-script scripts/compute_heatmap_data.rs out/hg002-k{k}-m{m}-kraken.bed out/heatmap-k{k}-m{m}-kraken.tsv")
        print(f"python3 scripts/draw_heatmap.py out/heatmap-k{k}-m{m}-kraken.tsv plots/heatmap-k{k}-m{m}-kraken.png --kraken")
