k = 63

for t in [1,2,4,8,16,32,64]:
    print(f"/usr/bin/time -v hks build -k 63 -o index/CHM13-k63.hks -t {t} --external-memory temp -f feature_tree/chromosome_fof.txt --color-names feature_tree/chromosome_names.txt --hierarchy feature_tree/chromosome_tree.txt -b CHM13/sbwt/CHM13-k63-added-dummies.sbwt -l CHM13/sbwt/CHM13-k63-added-dummies.lcs &> logs/hks-build-k{k}-t{t}.log")
