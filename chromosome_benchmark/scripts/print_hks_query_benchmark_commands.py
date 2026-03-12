for t in [1,2,4,8,16,32,64]:
    for k in [15,31,47,63]:
        print(f"/usr/bin/time -v hks lookup --report-misses --no-header -q query/hg002v1.1.fasta.gz -i index/CHM13-k63.hks -t {t} -k {k} 1> out/hg002-k{k}-t{t}-hks.tsv 2> logs/hks-query-k{k}-t{t}.log")

