for t in [1,2,4,8,16,32,64]:
    for k in [15,31,47,63]:
        print(f"/usr/bin/time -v hks lookup -q query/hg002v1.1_chr1.fasta -i index/CHM13-s63.hks -t {t} -k {k} 1> /dev/null 2> logs/hks-parallel-speedup-k{k}-t{t}.log")
