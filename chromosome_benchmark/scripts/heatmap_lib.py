"""Shared label-mapping utilities for draw_heatmap.py and draw_diff_heatmap.py."""

import numpy as np

HKS_HIERARCHY_FILE = "preprint_results/hks_hierarchy_dump.txt"
QUERY_NAMES_FILE   = "preprint_results/query_names.txt"
KRAKEN_NAMES_FILE  = "preprint_results/kraken_names.dmp"

# Fixed display order for classification labels (y-axis), top to bottom.
# invert_yaxis() is used, so ordered_rows is stored bottom-to-top (reversed).
ROW_ORDER_TOP_TO_BOTTOM = [
    "novel",
    "root", "autosome_multigroup1", "submetacentric_multigroup1",
    "metacentric_multigroup1", "acrocentric_multigroup1", "sex_multigroup1",
    "chrM", "chrY", "chrX",
    "chr22", "chr21", "chr20", "chr19", "chr18", "chr17", "chr16", "chr15",
    "chr14", "chr13", "chr12", "chr11", "chr10", "chr9", "chr8", "chr7",
    "chr6", "chr5", "chr4", "chr3", "chr2", "chr1",
]

# Fixed display order for query labels (x-axis), left to right.
COL_ORDER = [
    "chr1_MATERNAL", "chr1_PATERNAL", "chr2_MATERNAL", "chr2_PATERNAL",
    "chr3_MATERNAL", "chr3_PATERNAL", "chr4_MATERNAL", "chr4_PATERNAL",
    "chr5_MATERNAL", "chr5_PATERNAL", "chr6_MATERNAL", "chr6_PATERNAL",
    "chr7_MATERNAL", "chr7_PATERNAL", "chr8_MATERNAL", "chr8_PATERNAL",
    "chr9_MATERNAL", "chr9_PATERNAL", "chr10_MATERNAL", "chr10_PATERNAL",
    "chr11_MATERNAL", "chr11_PATERNAL", "chr12_MATERNAL", "chr12_PATERNAL",
    "chr13_MATERNAL", "chr13_PATERNAL", "chr14_MATERNAL", "chr14_PATERNAL",
    "chr15_MATERNAL", "chr15_PATERNAL", "chr16_MATERNAL", "chr16_PATERNAL",
    "chr17_MATERNAL", "chr17_PATERNAL", "chr18_MATERNAL", "chr18_PATERNAL",
    "chr19_MATERNAL", "chr19_PATERNAL", "chr20_MATERNAL", "chr20_PATERNAL",
    "chr21_MATERNAL", "chr21_PATERNAL", "chr22_MATERNAL", "chr22_PATERNAL",
    "chrX_MATERNAL", "chrY_PATERNAL", "chrM",
]


def load_hks_maps():
    """Return (class_taxid_to_name, query_taxid_to_name) for the HKS 0-indexed scheme.

    Reads HKS_HIERARCHY_FILE, which is a dump from --print-hierarchy: the first
    line is the number of nodes, followed by one label name per line in id order.
    """
    with open(HKS_HIERARCHY_FILE) as f:
        lines = [l.strip() for l in f if l.strip()]
    n = int(lines[0])
    class_taxid_to_name = {str(i): lines[1 + i] for i in range(n)}

    query_names = []
    with open(QUERY_NAMES_FILE) as f:
        for line in f:
            name = line.strip()
            if name:
                query_names.append(name)
    query_taxid_to_name = {str(i): name for i, name in enumerate(query_names)}

    return class_taxid_to_name, query_taxid_to_name


def load_kraken_taxid_map():
    taxid_to_name = {}
    with open(KRAKEN_NAMES_FILE) as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = parts[0].strip()
            name  = parts[1].strip()
            taxid_to_name[taxid] = name
    return taxid_to_name


def load_tsv(path):
    """Load a compute_heatmap_data TSV.

    Returns (tsv_query_labels, tsv_class_labels, mat) where mat has shape
    (n_queries, n_classes) with raw integer base counts.
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        tsv_class_labels = header[1:]
        tsv_query_labels = []
        data = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            tsv_query_labels.append(parts[0])
            data.append([int(x) for x in parts[1:]])
    return tsv_query_labels, tsv_class_labels, np.array(data, dtype=np.float64)


def translate(tsv_path, class_name_map, query_name_map, kraken_mode):
    """Load TSV, translate taxid labels to names, apply mode-specific fixes.

    In Kraken mode: taxid "0" -> "novel", "A" counts merged into "root".
    In HKS mode:    "-" -> "novel".

    Returns (col_labels, row_labels, mat) where:
      col_labels  -- query sequence names  (x-axis after transpose)
      row_labels  -- classification names  (y-axis after transpose)
      mat         -- shape (n_queries, n_classes), raw base counts
    """
    tsv_query_labels, tsv_class_labels, mat = load_tsv(tsv_path)

    col_labels = [query_name_map.get(l, l) for l in tsv_query_labels]
    row_labels = [class_name_map.get(l, l) for l in tsv_class_labels]

    if kraken_mode:
        row_labels = ["novel" if l in ("-", "0") else l for l in row_labels]
        if "A" in row_labels and "root" in row_labels:
            a_idx    = row_labels.index("A")
            root_idx = row_labels.index("root")
            mat[:, root_idx] += mat[:, a_idx]
            mat = np.delete(mat, a_idx, axis=1)
            row_labels.pop(a_idx)
    else:
        row_labels = ["novel" if l == "-" else l for l in row_labels]

    return col_labels, row_labels, mat
