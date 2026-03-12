for t in [64,32,16,8]:
    for k in [15,31,47,63]:
        m = 15 if k == 15 else 27
        p = 0 # For --minimizer-spaces.
        print(f"rm -f kraken_k{k}/hash.k2d") # If exists, Kraken does not build again
        print(f"/usr/bin/time -v kraken2-build --build --db kraken_k{k} --kmer-len {k} --threads {t} --minimizer-len {m} --minimizer-spaces {p} &> logs/kraken-build-k{k}-t{t}.log")
