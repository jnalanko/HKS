print("set -xue")
print("python3 scripts/setup_kraken_dir.py")
print("cp kraken_chromosome_db/taxonomy/names.dmp feature_tree/kraken_names.dmp")

k_values = [15,31,47,63]
m_values = [15,22,31] # Can be at most 31

for k in k_values:
    for m in m_values:
        if m <= k:
            dir = f"kraken_k{k}_m{m}"
            print(f"rm -rf {dir}")
            print(f"cp -r kraken_chromosome_db {dir}")

