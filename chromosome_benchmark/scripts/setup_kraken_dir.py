#!/usr/bin/env python3
"""Build Kraken2 taxonomy files and tagged FASTA library for human chromosomes."""

import os
import gzip

CHROM_NAMES_FILE = "feature_tree/chromosome_names.txt"
TREE_FILE = "feature_tree/chromosome_tree.txt"
CHROM_DIR = "CHM13/chromosomes"
OUT_DIR = "kraken_chromosome_db"

# --- Parse chromosome names (0-indexed: file i.fna = name at line i+1) ---
chrom_names = []
with open(CHROM_NAMES_FILE) as f:
    for line in f:
        name = line.strip()
        if name:
            chrom_names.append(name)

print(f"Loaded {len(chrom_names)} chromosome names: {chrom_names}")

# Taxid assignments: root=1, other internal nodes=2..n_internal, chromosomes=n_internal+1..
# (populated after internal nodes are parsed)

# --- Parse tree ---
with open(TREE_FILE) as f:
    lines = [l.strip() for l in f if l.strip()]

header = lines[0].split()
n_internal, n_edges = int(header[0]), int(header[1])

internal_nodes = []
for i in range(1, 1 + n_internal):
    internal_nodes.append(lines[i])

edges = []  # (child, parent)
for i in range(1 + n_internal, 1 + n_internal + n_edges):
    parts = lines[i].split()
    edges.append((parts[0], parts[1]))

print(f"Internal nodes: {internal_nodes}")
print(f"Edges: {edges[:5]} ...")

# Assign taxids: root=1, other internal nodes=2..n_internal, chromosomes=n_internal+1..
internal_taxid = {}
next_id = 2
for node in internal_nodes:
    if node == "root":
        internal_taxid[node] = 1
    else:
        internal_taxid[node] = next_id
        next_id += 1

chrom_taxid = {name: len(internal_nodes) + 1 + i for i, name in enumerate(chrom_names)}

def get_taxid(name):
    if name in internal_taxid:
        return internal_taxid[name]
    if name in chrom_taxid:
        return chrom_taxid[name]
    raise ValueError(f"Unknown node: {name}")

# --- Create output directories ---
os.makedirs(f"{OUT_DIR}/taxonomy", exist_ok=True)
os.makedirs(f"{OUT_DIR}/library/human", exist_ok=True)

# --- Write names.dmp ---
with open(f"{OUT_DIR}/taxonomy/names.dmp", "w") as f:
    def write_name(taxid, name):
        f.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")

    for node, taxid in internal_taxid.items():
        write_name(taxid, node)
    for name, taxid in chrom_taxid.items():
        write_name(taxid, name)

print(f"Wrote {OUT_DIR}/taxonomy/names.dmp")

# --- Write nodes.dmp ---
# Build parent map from edges
parent_map = {}
for child, parent in edges:
    parent_map[child] = parent
# root's parent is itself
parent_map["root"] = "root"

with open(f"{OUT_DIR}/taxonomy/nodes.dmp", "w") as f:
    def write_node(taxid, parent_taxid, rank):
        f.write(f"{taxid}\t|\t{parent_taxid}\t|\t{rank}\t|\n")

    # internal nodes (root first, parent of root is itself)
    for node, taxid in internal_taxid.items():
        parent_name = parent_map[node]
        write_node(taxid, get_taxid(parent_name), "no rank")
    # chromosomes
    for name, taxid in chrom_taxid.items():
        if name in parent_map:
            parent_name = parent_map[name]
            parent_tid = get_taxid(parent_name)
        else:
            # fallback: parent is root
            print(f"WARNING: {name} not in tree, parenting to root")
            parent_tid = internal_taxid["root"]
        rank = "mitochondrion" if name == "chrM" else "chromosome"
        write_node(taxid, parent_tid, rank)

print(f"Wrote {OUT_DIR}/taxonomy/nodes.dmp")

# --- Write tagged FASTA files and prelim_map.txt ---
prelim_map_path = os.path.join(OUT_DIR, "library", "human", "prelim_map.txt")
with open(prelim_map_path, "w") as prelim_map:
    for i, name in enumerate(chrom_names):
        src = os.path.join(CHROM_DIR, f"{i}.fna")
        dst = os.path.join(OUT_DIR, "library", "human", f"{name}.fna")
        taxid = chrom_taxid[name]

        if not os.path.exists(src):
            print(f"WARNING: {src} not found, skipping")
            continue

        print(f"Processing {src} -> {dst} (taxid={taxid})")
        with open(src, "rb") as fin, open(dst, "wb") as fout:
            for line in fin:
                if line.startswith(b">"):
                    # Insert kraken:taxid tag after '>'
                    rest = line[1:].rstrip()
                    fout.write(b">kraken:taxid|" + str(taxid).encode() + b"| " + rest + b"\n")
                    # Write seqid (first FASTA token) and taxid to prelim_map
                    seqid = f"kraken:taxid|{taxid}|"
                    prelim_map.write(f"TAXID\t{seqid}\t{taxid}\n")
                else:
                    fout.write(line)

print("Done.")
print(f"\nTaxid assignments:")
for node, taxid in sorted(internal_taxid.items(), key=lambda x: x[1]):
    print(f"  {taxid}: {node}")
for name, taxid in sorted(chrom_taxid.items(), key=lambda x: x[1]):
    print(f"  {taxid}: {name}")
