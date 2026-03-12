for t in [1,2,4,8,16,32,64]:
    for k in [15,31,47,63]:
        print(f"/usr/bin/time -v dks lookup -q query/hg002v1.1_chr1.fasta.gz -i index/CHM13-k63.dks -t {t} -k {k} 1> /dev/null 2> logs/dks-parallel-speedup-k{k}-t{t}.log")
