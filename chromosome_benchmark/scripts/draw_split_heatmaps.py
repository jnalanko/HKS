#!/usr/bin/env python3
"""Draw HKS and Kraken heatmaps in a 2×1 figure, and the diff heatmap as a separate figure.

Usage:
  python3 scripts/draw_split_heatmaps.py <hks.tsv> <kraken.tsv> <combined.png> <diff.png>
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

from heatmap_lib import (
    ROW_ORDER_TOP_TO_BOTTOM, COL_ORDER,
    load_hks_maps, load_kraken_taxid_map, translate,
)


def shorten_col_label(label):
    if label.endswith("_MATERNAL"):
        base = label[:-9]
        return (base[3:] if base.startswith("chr") else base) + "-M"
    if label.endswith("_PATERNAL"):
        base = label[:-9]
        return (base[3:] if base.startswith("chr") else base) + "-P"
    if label == "chrM":
        return "Mito"
    return label[3:] if label.startswith("chr") else label


def reorder(col_labels, row_labels, mat, ordered_cols, ordered_rows):
    col_idx = {l: i for i, l in enumerate(col_labels)}
    row_idx = {l: i for i, l in enumerate(row_labels)}
    q = [col_idx[l] for l in ordered_cols]
    r = [row_idx[l] for l in ordered_rows]
    return mat[np.ix_(q, r)]


def draw_log_panel(ax, fig, mat, ordered_cols, ordered_rows, title):
    log_mat = np.log10(mat + 1).T
    im = ax.imshow(log_mat, aspect="auto", cmap="inferno",
                   interpolation="nearest", origin="upper")
    ax.set_xticks(range(len(ordered_cols)))
    ax.set_xticklabels([shorten_col_label(l) for l in ordered_cols], rotation=90, ha="center", fontsize=32, fontfamily="monospace", fontweight="bold")
    ax.set_yticks(range(len(ordered_rows)))
    ax.set_yticklabels([l.replace("_multigroup1", "") for l in ordered_rows], fontsize=32, fontfamily="monospace", fontweight="bold")
    ax.invert_yaxis()
    ax.set_title(title, fontsize=80, pad=20)
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cbar.ax.tick_params(labelsize=32)
    cbar.set_label("log₁₀(count + 1)", fontsize=36)


def draw_diff_panel(ax, fig, mat, ordered_cols, ordered_rows):
    signed_log = np.sign(mat) * np.log10(np.abs(mat) + 1)
    im = ax.imshow(signed_log, aspect="auto", cmap="RdBu_r",
                   interpolation="nearest", origin="upper")
    ax.set_xticks(range(len(ordered_cols)))
    ax.set_xticklabels([shorten_col_label(l) for l in ordered_cols], rotation=90, ha="center", fontsize=32, fontfamily="monospace", fontweight="bold")
    ax.set_yticks(range(len(ordered_rows)))
    ax.set_yticklabels([l.replace("_multigroup1", "") for l in ordered_rows], fontsize=32, fontfamily="monospace", fontweight="bold")
    ax.invert_yaxis()
    ax.set_title("HKS minus Kraken", fontsize=80, pad=20)
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cbar.ax.tick_params(labelsize=32)
    cbar.set_label("sign(HKS−Kraken) · log₁₀(|HKS−Kraken| + 1)", fontsize=36)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hks_tsv",    help="HKS matrix TSV")
    parser.add_argument("kraken_tsv", help="Kraken matrix TSV")
    args = parser.parse_args()

    hks_class_map, query_name_map = load_hks_maps()
    kraken_class_map = load_kraken_taxid_map("feature_tree/kraken_names.dmp")

    hks_cols, hks_rows, hks_mat = translate(
        args.hks_tsv, hks_class_map, query_name_map, kraken_mode=False)
    kraken_cols, kraken_rows, kraken_mat = translate(
        args.kraken_tsv, kraken_class_map, query_name_map, kraken_mode=True)

    present_rows = set(hks_rows) | set(kraken_rows)
    ordered_rows = [l for l in reversed(ROW_ORDER_TOP_TO_BOTTOM) if l in present_rows]
    ordered_rows += sorted(present_rows - set(ROW_ORDER_TOP_TO_BOTTOM))

    present_cols = set(hks_cols) | set(kraken_cols)
    ordered_cols = [l for l in COL_ORDER if l in present_cols]
    ordered_cols += sorted(present_cols - set(COL_ORDER))

    hks_aligned    = reorder(hks_cols, hks_rows, hks_mat, ordered_cols, ordered_rows)
    kraken_aligned = reorder(kraken_cols, kraken_rows, kraken_mat, ordered_cols, ordered_rows)
    diff = (hks_aligned - kraken_aligned).T

    w = max(8, len(ordered_cols) * 0.6)
    h = max(6, len(ordered_rows) * 0.35)

    # Figure 1: HKS and Kraken side by side (2 columns, 1 row)
    fig1, axes1 = plt.subplots(1, 2, figsize=(w * 2, h * 2))
    draw_log_panel(axes1[0], fig1, hks_aligned,    ordered_cols, ordered_rows, "HKS")
    draw_log_panel(axes1[1], fig1, kraken_aligned, ordered_cols, ordered_rows, "Kraken")
    fig1.tight_layout()
    fig1.savefig("plots/heatmaps.pdf")
    print("Saved plots/heatmaps.pdf")

    # Figure 2: diff heatmap alone
    fig2, ax2 = plt.subplots(1, 1, figsize=(w, h * 2))
    draw_diff_panel(ax2, fig2, diff, ordered_cols, ordered_rows)
    fig2.tight_layout()
    fig2.savefig("plots/diff_heatmap.pdf")
    print("Saved plots/diff_heatmap.pdf")


if __name__ == "__main__":
    main()
