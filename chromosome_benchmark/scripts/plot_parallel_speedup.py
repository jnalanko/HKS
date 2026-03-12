#!/usr/bin/env python3
"""Parse parallel speedup logs and plot throughput + speedup."""

from __future__ import annotations
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

LOG_DIR = Path("final_logs/parallel_speedup")
K_VALUES = [15, 31, 47, 63]
#T_VALUES = [1, 2, 4, 8, 16, 32, 64]
T_VALUES = [1, 2, 4, 8, 16, 32]

TIME_RE = re.compile(r"Query time per base pair:\s+([\d.]+)\s+nanoseconds")


def parse_ns_per_bp(path: Path) -> float:
    text = path.read_text()
    m = TIME_RE.search(text)
    if m is None:
        raise ValueError(f"No timing line found in {path}")
    return float(m.group(1))


def load_data() -> dict[int, dict[int, float]]:
    """Returns {k: {t: ns_per_bp}}"""
    data: dict[int, dict[int, float]] = {}
    missing = []
    for k in K_VALUES:
        data[k] = {}
        for t in T_VALUES:
            path = LOG_DIR / f"dks-parallel-speedup-k{k}-t{t}.log"
            if not path.exists():
                missing.append(str(path))
                continue
            data[k][t] = parse_ns_per_bp(path)
    if missing:
        print(f"Warning: {len(missing)} log file(s) not found:", file=sys.stderr)
        for p in missing:
            print(f"  {p}", file=sys.stderr)
    return data


def main():
    data = load_data()

    fig, (ax_thr, ax_spd) = plt.subplots(1, 2, figsize=(12, 5))
    colors = plt.cm.tab10(np.linspace(0, 0.4, len(K_VALUES)))

    for color, k in zip(colors, K_VALUES):
        kdata = data[k]
        if not kdata:
            continue
        ts = sorted(kdata)
        ns = [kdata[t] for t in ts]
        throughput = [1e3 / n for n in ns]  # Mbases/s

        ax_thr.plot(ts, throughput, marker="o", label=f"s = {k}, k = 63", color=color)

        # Speedup relative to single-threaded
        if 1 in kdata:
            baseline = kdata[1]
            speedup = [baseline / kdata[t] for t in ts]
            ax_spd.plot(ts, speedup, marker="o", label=f"s = {k}, k = 63", color=color)

    # Ideal speedup line
    ax_spd.plot(T_VALUES, T_VALUES, "k--", linewidth=1, label="ideal")

    # --- left panel ---
    ax_thr.set_xlabel("Threads (t)")
    ax_thr.set_ylabel("Throughput (Mbases / s)")
    ax_thr.set_title("Throughput vs threads")
    ax_thr.set_xticks(T_VALUES)
    ax_thr.set_ylim(bottom=0)
    ax_thr.legend()
    ax_thr.grid(True, alpha=0.3)

    # --- right panel ---
    ax_spd.set_xlabel("Threads (t)")
    ax_spd.set_ylabel("Speedup vs t=1")
    ax_spd.set_title("Parallel speedup")
    ax_spd.set_xticks(T_VALUES)
    ax_spd.set_ylim(bottom=0)
    ax_spd.legend()
    ax_spd.grid(True, alpha=0.3)

    fig.tight_layout()
    out = Path("plots/parallel_speedup.pdf")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
