#!/usr/bin/env python3
"""Benchmark POST /api/run latency and server memory vs tree size.

Usage:
    python outputs/bench/benchmark_api.py [--port 8765] [--repeats 5]

Outputs:
    outputs/bench/bench_results.csv
"""

import argparse
import csv
import io
import os
import statistics
import subprocess
import sys
import time

import random

import dendropy
import psutil
import requests

# Tree sizes to benchmark (number of taxa)
TREE_SIZES = [100, 500, 1000, 2000, 5000]
SEED = 42
N_REPS = 5        # requests per tree size
FIXED_N = 20      # representatives to select
WARMUP = 2        # warm-up requests before timing

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
OUT_CSV = os.path.join(SCRIPT_DIR, "bench_results.csv")


# ---------------------------------------------------------------------------
# Tree generation
# ---------------------------------------------------------------------------

def generate_newick(n_taxa: int, seed: int) -> str:
    """Return a Newick string for a birth-death tree with n_taxa leaves."""
    rng = random.Random(seed)
    tree = dendropy.simulate.treesim.birth_death_tree(
        birth_rate=1.0,
        death_rate=0.0,
        num_extant_tips=n_taxa,
        rng=rng,
    )
    return tree.as_string(schema="newick")


# ---------------------------------------------------------------------------
# Server management
# ---------------------------------------------------------------------------

def start_server(port: int) -> subprocess.Popen:
    """Start parnas server as a subprocess via uv run."""
    env = os.environ.copy()
    env.pop("VIRTUAL_ENV", None)  # don't let parent venv override uv's choice
    proc = subprocess.Popen(
        ["uv", "run", "parnas", "server", "--port", str(port)],
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
    )
    # Wait until ready (poll /api/run endpoint)
    url = f"http://localhost:{port}"
    # Wait for server to start (only break on successful connection from OUR process)
    for _ in range(90):
        if proc.poll() is not None:
            out, _ = proc.communicate()
            raise RuntimeError(f"Server exited early: {out[:300]}")
        try:
            requests.get(url, timeout=1)
            break
        except Exception:
            time.sleep(0.5)
    else:
        proc.terminate()
        raise RuntimeError(f"Server did not start on port {port}")
    time.sleep(1)  # let warmup solver finish
    return proc


def stop_server(proc: subprocess.Popen) -> None:
    proc.terminate()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()


# ---------------------------------------------------------------------------
# Single request benchmark
# ---------------------------------------------------------------------------

def post_run(url: str, newick: str) -> float:
    """POST /api/run and return wall-clock seconds."""
    data = {
        "cover": "false",
        "evaluate": "false",
        "binary": "false",
        "sweep": "none",
        "n": str(FIXED_N),
    }
    files = {
        "tree": ("tree.nwk", io.BytesIO(newick.encode()), "text/plain"),
    }
    t0 = time.perf_counter()
    resp = requests.post(f"{url}/api/run", data=data, files=files, timeout=300)
    elapsed = time.perf_counter() - t0
    if resp.status_code != 200:
        raise RuntimeError(f"HTTP {resp.status_code}: {resp.text[:200]}")
    return elapsed


def sample_server_rss(proc: subprocess.Popen) -> float:
    """Return current RSS in MB for the server process (best-effort)."""
    try:
        ps = psutil.Process(proc.pid)
        return ps.memory_info().rss / 1024 / 1024
    except Exception:
        return float("nan")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--repeats", type=int, default=N_REPS)
    parser.add_argument("--sizes", nargs="+", type=int, default=TREE_SIZES)
    args = parser.parse_args()

    url = f"http://localhost:{args.port}"

    print(f"Starting PARNAS server on port {args.port} ...", flush=True)
    proc = start_server(args.port)
    print("Server ready.", flush=True)

    rows = []
    try:
        for n_taxa in args.sizes:
            print(f"\n--- n_taxa={n_taxa} ---", flush=True)
            newick = generate_newick(n_taxa, SEED)

            # Warm-up
            for _ in range(WARMUP):
                post_run(url, newick)

            latencies = []
            rss_samples = []
            for rep in range(args.repeats):
                rss_before = sample_server_rss(proc)
                lat = post_run(url, newick)
                rss_after = sample_server_rss(proc)
                peak_mb = max(rss_before, rss_after)
                latencies.append(lat)
                rss_samples.append(peak_mb)
                print(
                    f"  rep {rep+1}: latency={lat:.3f}s  RSS≈{peak_mb:.1f}MB",
                    flush=True,
                )

            median_lat = statistics.median(latencies)
            iqr_lat = (
                sorted(latencies)[int(len(latencies) * 0.75)]
                - sorted(latencies)[int(len(latencies) * 0.25)]
            )
            peak_mem = max(rss_samples)
            rows.append(
                {
                    "n_taxa": n_taxa,
                    "latency_median_s": round(median_lat, 4),
                    "latency_iqr_s": round(iqr_lat, 4),
                    "peak_mem_mb": round(peak_mem, 1),
                }
            )
            print(
                f"  → median={median_lat:.3f}s  IQR={iqr_lat:.3f}s  peakMem={peak_mem:.1f}MB"
            )
    finally:
        stop_server(proc)
        print("\nServer stopped.", flush=True)

    # Write CSV
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["n_taxa", "latency_median_s", "latency_iqr_s", "peak_mem_mb"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nResults written to {OUT_CSV}")
    for row in rows:
        print(
            f"  n_taxa={row['n_taxa']:>6}  lat={row['latency_median_s']}s"
            f"  mem={row['peak_mem_mb']}MB"
        )


if __name__ == "__main__":
    main()
