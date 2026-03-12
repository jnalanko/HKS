#!/usr/bin/env python3
"""Draw a heatmap from a compute_heatmap_data TSV with logarithmic brightness.

Label order is derived from feature_tree/chromosome_names.txt and
feature_tree/chromosome_tree.txt. Labels not found in those files are
appended at the end in lexicographic order.

Usage:
  python3 draw_heatmap.py <matrix.tsv> <out.png>
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

from heatmap_lib import (
    ROW_ORDER_TOP_TO_BOTTOM, COL_ORDER,
    load_hks_maps, load_kraken_taxid_map, translate,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tsv", help="matrix TSV from compute_heatmap_data")
    parser.add_argument("out", help="output PNG path")
    parser.add_argument("--kraken", action="store_true",
                        help="resolve taxids using feature_tree/kraken_names.dmp")
    args = parser.parse_args()

    hks_class_map, query_taxid_to_name = load_hks_maps()
    class_name_map = load_kraken_taxid_map("feature_tree/kraken_names.dmp") if args.kraken else hks_class_map

    col_labels, row_labels, mat = translate(args.tsv, class_name_map, query_taxid_to_name, args.kraken)

    # Determine display order from fixed label lists; append any unlisted labels at the end.
    present_rows = set(row_labels)
    ordered_rows = [l for l in reversed(ROW_ORDER_TOP_TO_BOTTOM) if l in present_rows]
    ordered_rows += sorted(present_rows - set(ROW_ORDER_TOP_TO_BOTTOM))

    present_cols = set(col_labels)
    ordered_cols = [l for l in COL_ORDER if l in present_cols]
    ordered_cols += sorted(present_cols - set(COL_ORDER))

    # Build index maps into the TSV orientation (mat rows=queries, mat cols=classifications)
    query_idx = {l: i for i, l in enumerate(col_labels)}
    class_idx = {l: i for i, l in enumerate(row_labels)}

    # Reorder: select query rows in ordered_cols order, class cols in ordered_rows order
    q_indices = [query_idx[l] for l in ordered_cols]
    c_indices = [class_idx[l] for l in ordered_rows]
    mat = mat[np.ix_(q_indices, c_indices)]

    log_mat = np.log10(mat + 1).T  # transpose → (n_ordered_rows, n_ordered_cols)

    fig, ax = plt.subplots(figsize=(max(8, len(ordered_cols) * 0.5),
                                    max(6, len(ordered_rows) * 0.35)))
    im = ax.imshow(log_mat, aspect="auto", cmap="inferno",
                   interpolation="nearest", origin="upper")

    ax.set_xticks(range(len(ordered_cols)))
    ax.set_xticklabels(ordered_cols, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(ordered_rows)))
    ax.set_yticklabels(ordered_rows, fontsize=7)
    ax.invert_yaxis()

    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cbar.set_label("log₁₀(bases + 1)", fontsize=9)

    ax.set_xlabel("query label")
    ax.set_ylabel("classification label")

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved {args.out}")

    print("\nRow labels (y-axis, top to bottom):")
    for i, label in enumerate(reversed(ordered_rows)):
        print(f"  {i}: {label}")
    print("\nColumn labels (x-axis, left to right):")
    for i, label in enumerate(ordered_cols):
        print(f"  {i}: {label}")


if __name__ == "__main__":
    main()
