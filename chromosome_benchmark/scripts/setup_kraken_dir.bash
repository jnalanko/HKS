set -xue

python3 scripts/setup_kraken_dir.py

rm -rf kraken_k15
rm -rf kraken_k31
rm -rf kraken_k47
rm -rf kraken_k63

cp -r kraken_chromosome_db kraken_k15
cp -r kraken_chromosome_db kraken_k31
cp -r kraken_chromosome_db kraken_k47
cp -r kraken_chromosome_db kraken_k63
cp kraken_chromosome_db/taxonomy/names.dmp feature_tree/kraken_names.dmp
