s = 63

for t in [1,2,4,8,16,32,64]:
    print(f"/usr/bin/time -v hks build -s 63 -o index/CHM13-s63.hks -t {t} --external-memory temp --label-by-file feature_tree/chromosome_fof.txt --labels feature_tree/chromosome_names.txt --hierarchy feature_tree/chromosome_tree.txt --load-sbwt CHM13/sbwt/CHM13-k63-added-dummies.sbwt --load-lcs CHM13/sbwt/CHM13-k63-added-dummies.lcs &> logs/hks-build-s{s}-t{t}.log")


#--label-by-file feature_tree/chromosome_fof.txt --labels feature_tree/chromosome_names.txt
