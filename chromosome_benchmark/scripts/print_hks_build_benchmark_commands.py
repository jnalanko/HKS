s = 63

for t in [1,2,4,8,16,32,64]:
    build_base = (
        f"hks build-base -s {s} -o index/CHM13-s{s}.hksb -t {t} "
        f"--external-memory temp "
        f"--load-sbwt CHM13/sbwt/CHM13-k63-added-dummies.sbwt "
        f"--load-lcs CHM13/sbwt/CHM13-k63-added-dummies.lcs"
    )
    add_feature_set = (
        f"hks add-feature-set -i index/CHM13-s{s}.hksb -o index/CHM13-s{s}.hksf "
        f"--feature-file-list feature_tree/chromosome_fof.txt "
        f"--feature-names feature_tree/chromosome_names.txt "
        f"--feature-hierarchy feature_tree/chromosome_tree.txt "
        f"--feature-set-name chromosome --variable-k-support -t {t}"
    )
    print(f'/usr/bin/time -v bash -c "{build_base} && {add_feature_set}" &> logs/hks-build-s{s}-t{t}.log')
