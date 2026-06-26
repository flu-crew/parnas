#!/usr/bin/env python3
"""Generate benchmark figure and LaTeX table from bench_results.csv.

Usage:
    python outputs/bench/plot_bench.py

Outputs:
    manuscript/images/bench-latency.png
    outputs/bench/bench_table.tex
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

CSV_PATH = os.path.join(SCRIPT_DIR, "bench_results.csv")
IMG_OUT = os.path.join(SCRIPT_DIR, "bench-latency.png")

def load_results(path: str):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "n_taxa": int(row["n_taxa"]),
                "latency_median_s": float(row["latency_median_s"]),
                "latency_iqr_s": float(row["latency_iqr_s"]),
                "peak_mem_mb": float(row["peak_mem_mb"]),
            })
    return rows


def make_plot(rows, out_path: str):
    n = [r["n_taxa"] for r in rows]
    lat = [r["latency_median_s"] for r in rows]
    iqr = [r["latency_iqr_s"] for r in rows]
    mem = [r["peak_mem_mb"] for r in rows]

    fig, ax1 = plt.subplots(figsize=(5, 3.5))

    color_lat = "#2166ac"
    color_mem = "#d6604d"

    ax1.errorbar(n, lat, yerr=[i / 2 for i in iqr],
                 color=color_lat, marker="o", linewidth=1.5,
                 capsize=3, label="API latency (s)")
    ax1.set_xlabel("Number of taxa")
    ax1.set_ylabel("Median API latency (s)", color=color_lat)
    ax1.tick_params(axis="y", labelcolor=color_lat)
    ax1.set_xscale("log")
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax1.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)

    ax2 = ax1.twinx()
    ax2.plot(n, mem, color=color_mem, marker="s", linestyle="--",
             linewidth=1.5, label="Peak memory (MB)")
    ax2.set_ylabel("Peak server RSS (MB)", color=color_mem)
    ax2.tick_params(axis="y", labelcolor=color_mem)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc="upper left")

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Figure written to {out_path}")

def main():
    if not os.path.exists(CSV_PATH):
        print(f"ERROR: {CSV_PATH} not found. Run benchmark_api.py first.")
        raise SystemExit(1)
    rows = load_results(CSV_PATH)
    make_plot(rows, IMG_OUT)


if __name__ == "__main__":
    main()
