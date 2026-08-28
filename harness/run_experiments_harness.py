#!/usr/bin/env python3
"""
run_experiments_harness.py

Generates synthetic FASTA sequence pairs, runs two alignment binaries
(tiled_sw.exe and wfa_banded.exe) across parameter sweeps, collects timings &
memory use, writes CSV files, and produces PNG plots.

Usage example (from harness folder):
python run_experiments_harness.py --outdir ..\results --pairs-per-identity 200 \
    --lengths 2000 5000 --identities 0.98 0.95 0.90 0.80 --threads-list 1 2 4 8 \
    --repeats 3 --auto-compile
"""

import argparse
import csv
import os
import random
import subprocess
import sys
import time
from pathlib import Path
from statistics import median, mean, stdev

# Try optional dependencies
try:
    import psutil
except:
    psutil = None

try:
    import matplotlib.pyplot as plt
    import numpy as np
except:
    plt = None
    np = None


# ============================================================
# Utility: Mutation (identity-based)
# ============================================================
def mutate(seq, identity):
    out = []
    keep = int(len(seq) * identity)
    idx = set(random.sample(range(len(seq)), keep))

    for i, c in enumerate(seq):
        if i in idx:
            out.append(c)
        else:
            out.append(random.choice([x for x in "ACGT" if x != c]))

    return "".join(out)


# ============================================================
# FASTA writer
# ============================================================
def write_fasta_pairs(path: Path, pairs):
    with path.open("w", newline="") as fh:
        for (ida, sa, idb, sb) in pairs:
            fh.write(f">{ida}\n{sa}\n>{idb}\n{sb}\n")


# ============================================================
# Run binary and collect metrics
# ============================================================
def call_binary(cmd, timeout=None):
    start = time.perf_counter()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    peak_kb = None
    if psutil:
        p = psutil.Process(proc.pid)
    else:
        p = None

    # Poll process
    try:
        while True:
            ret = proc.poll()
            if ret is not None:
                break
            if p:
                try:
                    rss = p.memory_info().rss
                    if peak_kb is None or rss > peak_kb:
                        peak_kb = rss
                except:
                    pass
            time.sleep(0.01)

        # last memory check
        if p:
            try:
                rss = p.memory_info().rss
                if peak_kb is None or rss > peak_kb:
                    peak_kb = rss
            except:
                pass

        out, err = proc.communicate(timeout=1)

    except subprocess.TimeoutExpired:
        proc.kill()
        out, err = proc.communicate()

    end = time.perf_counter()

    # convert memory
    if peak_kb is not None:
        peak_kb = peak_kb // 1024

    # decode
    out = out.decode("utf-8", errors="replace")
    err = err.decode("utf-8", errors="replace")

    return end - start, peak_kb, out.strip(), err.strip()


# ============================================================
# Generate pairs
# ============================================================
def generate_pairs_for_identity(n, length, identity):
    base = ''.join(random.choice("ACGT") for _ in range(length))
    pairs = []
    for i in range(n):
        ida = f"pair{i}_a"
        idb = f"pair{i}_b"
        sa = base
        sb = mutate(sa, identity)
        pairs.append((ida, sa, idb, sb))
    return pairs


# ============================================================
# Main sweep
# ============================================================
def run_sweep(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw_csv = outdir / "raw_runs.csv"
    agg_csv = outdir / "aggregated_results.csv"

    fieldnames = [
        "binary","length","identity","threads","repeat",
        "pairs_file","n_pairs","elapsed_s","peak_rss_kb",
        "stdout","stderr"
    ]

    raw_rows = []

    with raw_csv.open("w", newline="", encoding="utf-8") as rf:
        w = csv.DictWriter(rf, fieldnames=fieldnames)
        w.writeheader()

        for L in args.lengths:
            for ident in args.identities:
                print(f"[info] Generating pairs: length={L}, identity={ident}")
                pairs = generate_pairs_for_identity(args.pairs_per_identity, L, ident)
                fasta_path = outdir / f"pairs_L{L}_id{int(ident*100)}.fa"
                write_fasta_pairs(fasta_path, pairs)

                for t in args.threads_list:
                    for r in range(args.repeats):

                        for binary in args.binaries:
                            bin_path = Path(binary)
                            if not bin_path.is_absolute():
                                bin_path = (Path(args.root_dir) / binary).resolve()

                            if not bin_path.exists():
                                print(f"[ERROR] Binary not found: {bin_path}")
                                continue

                            # build command
                            if "tiled" in bin_path.name.lower():
                                cmd = [str(bin_path), str(fasta_path), str(args.tile_size), str(t)]
                            else:
                                cmd = [str(bin_path), str(fasta_path), str(args.band_threshold), str(t)]

                            print("[run]", cmd)

                            elapsed, peak, out, err = call_binary(cmd, timeout=args.timeout)

                            row = {
                                "binary": bin_path.name,
                                "length": L,
                                "identity": ident,
                                "threads": t,
                                "repeat": r,
                                "pairs_file": str(fasta_path),
                                "n_pairs": args.pairs_per_identity,
                                "elapsed_s": elapsed,
                                "peak_rss_kb": peak if peak else "",
                                "stdout": out[:6000],
                                "stderr": err[:2000],
                            }
                            w.writerow(row)
                            raw_rows.append(row)

    # -----------------------------
    # Aggregate
    # -----------------------------
    groups = {}
    for row in raw_rows:
        key = (row["binary"], row["length"], row["identity"], row["threads"])
        groups.setdefault(key, []).append(row)

    with agg_csv.open("w", newline="", encoding="utf-8") as af:
        fn2 = ["binary","length","identity","threads","n_runs",
               "median_time","mean_time","stdev_time","median_rss_kb"]
        w2 = csv.DictWriter(af, fieldnames=fn2)
        w2.writeheader()

        for key, rows in sorted(groups.items()):
            times = [r["elapsed_s"] for r in rows]
            rss = [int(r["peak_rss_kb"]) for r in rows if r["peak_rss_kb"]]
            rec = {
                "binary": key[0],
                "length": key[1],
                "identity": key[2],
                "threads": key[3],
                "n_runs": len(rows),
                "median_time": median(times),
                "mean_time": mean(times),
                "stdev_time": stdev(times) if len(times) > 1 else 0.0,
                "median_rss_kb": median(rss) if rss else "",
            }
            w2.writerow(rec)

    print("[done] Wrote:", raw_csv)
    print("[done] Wrote:", agg_csv)

    # -----------------------------
    # Plots
    # -----------------------------
    if plt is None:
        print("[warn] matplotlib not installed — skipping plots")
        return

    # Load aggregated into pandas
    import pandas as pd
    df = pd.read_csv(agg_csv)

    # Runtime vs identity (threads=1)
    fig = plt.figure(figsize=(7,4))
    for bname, g in df[df.threads==1].groupby("binary"):
        g2 = g.sort_values("identity", ascending=False)
        plt.plot(g2.identity*100, g2.median_time, marker='o', label=bname)
    plt.gca().invert_xaxis()
    plt.xlabel("Identity (%)")
    plt.ylabel("Median Runtime (s)")
    plt.title("Runtime vs Identity (threads=1)")
    plt.grid(True)
    plt.legend()
    fig.tight_layout()
    p1 = outdir / "runtime_vs_identity.png"
    fig.savefig(p1, dpi=200)
    plt.close(fig)
    print("[done] Plot:", p1)

    # Scaling plot for largest length
    Lmax = df.length.max()
    for ident in sorted(df.identity.unique(), reverse=True):
        fig = plt.figure(figsize=(7,4))
        tmp = df[(df.length==Lmax) & (df.identity==ident)]
        for bname, g in tmp.groupby("binary"):
            g2 = g.sort_values("threads")
            plt.plot(g2.threads, g2.median_time, marker='o', label=bname)
        plt.xlabel("Threads")
        plt.ylabel("Median Runtime (s)")
        plt.title(f"Scaling at length={Lmax}, identity={int(ident*100)}%")
        plt.grid(True)
        plt.legend()
        fig.tight_layout()
        p2 = outdir / f"scaling_L{Lmax}_id{int(ident*100)}.png"
        fig.savefig(p2, dpi=200)
        plt.close(fig)
        print("[done] Plot:", p2)


# ============================================================
# Auto-compile (optional)
# ============================================================
def auto_compile(args):
    srcdir = Path(args.root_dir) / "src"
    tiled = srcdir / "tiled_sw.c"
    wfa1 = srcdir / "wfa_bounded.c"
    wfa2 = srcdir / "wfa_banded.c"

    print("[compile] Attempting auto-compile…")

    if tiled.exists():
        print("[compile] tiled_sw.c")
        subprocess.call(["gcc","-O3","-fopenmp",str(tiled),"-o",str(Path(args.root_dir)/"tiled_sw.exe")])

    if wfa1.exists():
        print("[compile] wfa_bounded.c")
        subprocess.call(["gcc","-O3","-fopenmp",str(wfa1),"-o",str(Path(args.root_dir)/"wfa_banded.exe")])
    elif wfa2.exists():
        print("[compile] wfa_banded.c")
        subprocess.call(["gcc","-O3","-fopenmp",str(wfa2),"-o",str(Path(args.root_dir)/"wfa_banded.exe")])



# ============================================================
# CLI
# ============================================================
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="../results")
    p.add_argument("--pairs-per-identity", type=int, default=200)
    p.add_argument("--lengths", nargs="+", type=int, default=[2000,5000])
    p.add_argument("--identities", nargs="+", type=float, default=[0.98,0.95,0.90,0.80])
    p.add_argument("--threads-list", nargs="+", type=int, default=[1,2,4,8])
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--tile-size", type=int, default=128)
    p.add_argument("--band-threshold", type=int, default=200)
    p.add_argument("--binaries", nargs="+",
                  default=["..\\tiled_sw.exe", "..\\wfa_banded.exe"])
    p.add_argument("--root-dir", default="..")
    p.add_argument("--auto-compile", action="store_true")

    args = p.parse_args()

    if args.auto_compile:
        auto_compile(args)

    run_sweep(args)


if __name__ == "__main__":
    main()
