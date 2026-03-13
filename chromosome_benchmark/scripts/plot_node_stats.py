#!/usr/bin/env python3
"""Plot node_stats.tsv: x=s, y=count, one line per color.

Usage:
  python3 plot_node_stats.py
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

TSV = "preprint_results/node_stats.tsv"
OUT = "plots/node_stats.pdf"

plt.rcParams.update({
    "axes.labelsize":  20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
})


INTERNAL_ORDER = [
    "root",
    "autosome_multigroup1",
    "sex_multigroup1",
    "acrocentric_multigroup1",
    "metacentric_multigroup1",
    "submetacentric_multigroup1",
]

def node_sort_key(name):
    if name in INTERNAL_ORDER:
        return (0, INTERNAL_ORDER.index(name), 0)
    if name.startswith("chr"):
        suffix = name[3:]
        if suffix.isdigit():
            return (1, int(suffix), 0)
        return (2, {"X": 0, "Y": 1, "M": 2}.get(suffix, 99), 0)
    return (3, 0, 0)

df = pd.read_csv(TSV, sep="\t")

all_nodes = sorted(df["color"].unique(), key=node_sort_key)
palette = [cm.tab20(i) for i in range(20)] + [cm.tab20b(i) for i in range(20)]
node_color = {node: palette[i] for i, node in enumerate(all_nodes)}

total_per_s = df.groupby("s")["count"].sum().rename("total")
df = df.join(total_per_s, on="s")
df["pct"] = df["count"] / df["total"] * 100

fig, axes = plt.subplots(1, 2, figsize=(18, 6))

for color in all_nodes:
    group = df[df["color"] == color].sort_values("s")
    axes[0].plot(group["s"].to_numpy(), group["count"].to_numpy(),
                 color=node_color[color], label=color.replace("_multigroup1", ""))

s_values = sorted(df["s"].unique())
pct_matrix = []
labels = []
colors = []
for color in all_nodes:
    group = df[df["color"] == color].set_index("s").reindex(s_values)["pct"].fillna(0)
    pct_matrix.append(group.to_numpy())
    labels.append(color.replace("_multigroup1", ""))
    colors.append(node_color[color])

axes[1].stackplot(s_values, pct_matrix, labels=labels, colors=colors)

axes[0].set_xlabel(r"$k$")
axes[0].set_ylabel("Count")
axes[0].set_xlim(1, 63)
axes[1].set_xlabel(r"$k$")
axes[1].set_ylabel("% of Total")
axes[1].set_xlim(1, 63)
axes[1].set_ylim(0, 100)
fig.subplots_adjust(bottom=0.35)
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.015),
           ncol=10, fontsize=14)
plt.savefig(OUT, dpi=150)
print(f"Saved {OUT}")
