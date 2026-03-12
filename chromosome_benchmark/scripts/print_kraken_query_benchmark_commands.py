for t in [64,32,16,8,4,2,1]:
    for k in [15,31,47,63]:
        print(f"/usr/bin/time -v kraken2 --db kraken_k{k} --gzip-compressed --output out/hg002-k{k}-t{t}-kraken.txt --threads {t} query/hg002v1.1.fasta.gz &> logs/kraken-query-k{k}-t{t}.log")
