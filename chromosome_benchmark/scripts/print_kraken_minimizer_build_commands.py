k_values = [15,31,47,63]
m_values = [15,22,31] # Can be at most 31
t = 32

print("set -x")
for k in k_values:
    for m in m_values:
        if m <= k:
            p = 0 # For --minimizer-spaces.
            print(f"rm -f kraken_k{k}_m{m}/hash.k2d") # If exists, Kraken does not build again
            print(f"/usr/bin/time -v kraken2-build --build --db kraken_k{k}_m{m} --kmer-len {k} --threads {t} --minimizer-len {m} --minimizer-spaces {p} &> logs/kraken-build-k{k}-m{m}-t{t}.log")
