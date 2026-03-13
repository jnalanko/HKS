k_values = [15,31,47,63]
m_values = [15,22,31] # Can be at most 31
t = 32

print("set -x")


for k in k_values:
    for m in m_values:
        print(f"/usr/bin/time -v kraken2 --db kraken_k{k}_m{m} --gzip-compressed --output out/hg002-k{k}-m{m}-t{t}-kraken.txt --threads {t} query/hg002v1.1.fasta.gz &> logs/kraken-query-k{k}-m{m}-t{t}.log")

