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
import ctypes
import csv
import os
import re
import random
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from statistics import median, mean, stdev

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PERF_EVENTS = "L1-dcache-load-misses,instructions,cycles"

# Try optional dependencies
try:
    import psutil
except ImportError:
    psutil = None

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


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
def parse_perf_metrics(stderr):
    metrics = {"l1_dcache_load_misses": "", "instructions": "", "cycles": ""}
    event_names = {
        "l1-dcache-load-misses": "l1_dcache_load_misses",
        "instructions": "instructions",
        "cycles": "cycles",
    }

    for line in stderr.splitlines():
        fields = line.strip().split(",")
        if len(fields) < 3:
            continue
        event = fields[2].strip().lower()
        metric = event_names.get(event)
        if metric and re.fullmatch(r"\d+", fields[0].strip()):
            metrics[metric] = int(fields[0].strip())

    if metrics["instructions"] != "" and metrics["cycles"] not in ("", 0):
        metrics["ipc"] = metrics["instructions"] / metrics["cycles"]
    else:
        metrics["ipc"] = ""
    return metrics


def set_windows_process_affinity(process, affinity_mask):
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    process_handle = ctypes.c_void_p(process._handle)

    set_affinity = kernel32.SetProcessAffinityMask
    set_affinity.argtypes = (ctypes.c_void_p, ctypes.c_size_t)
    set_affinity.restype = ctypes.c_int
    if not set_affinity(process_handle, ctypes.c_size_t(affinity_mask)):
        raise ctypes.WinError(ctypes.get_last_error())

def parse_vtune_report(report_path):
    metrics = {"l1_dcache_load_misses": "", "instructions": "", "cycles": "", "ipc": ""}
    labels = {
        "l1_dcache_load_misses": (
            "l1 data cache load miss", "l1-dcache-load-miss",
            "l1d cache miss", "l1 data cache miss",
        ),
        "instructions": (
            "instructions retired", "retired instructions", "instructions executed",
        ),
        "cycles": ("cpu cycles", "unhalted core cycles", "clockticks"),
        "ipc": ("instructions per cycle", "instructions per clock", "ipc"),
    }

    def numeric_value(value):
        cleaned = value.strip().replace(",", "")
        match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", cleaned)
        return float(match.group(0)) if match else None

    def row_value(row):
        if ":" in row[0]:
            value = row[0].split(":", 1)[1]
            if len(row) > 1:
                value += "," + ",".join(row[1:])
            return numeric_value(value)
        for index, value in enumerate(row[1:], start=1):
            parsed = numeric_value(value)
            if parsed is None:
                continue
            if parsed.is_integer():
                digits = value.strip()
                combined = digits
                next_index = index + 1
                while (
                    next_index < len(row)
                    and re.fullmatch(r"\d{3}", row[next_index].strip())
                ):
                    combined += row[next_index].strip()
                    next_index += 1
                if combined != digits:
                    parsed = float(combined)
            return parsed
        return None

    with Path(report_path).open(newline="", encoding="utf-8-sig", errors="replace") as report:
        for row in csv.reader(report):
            if not row:
                continue
            normalized = " ".join(row).lower()
            for metric, names in labels.items():
                if metrics[metric] != "":
                    continue
                if any(name in normalized for name in names):
                    parsed = row_value(row)
                    if parsed is not None:
                        metrics[metric] = parsed

    if any(metrics[key] == "" for key in metrics):
        text = Path(report_path).read_text(encoding="utf-8-sig", errors="replace")
        for metric, names in labels.items():
            if metrics[metric] != "":
                continue
            for line in text.splitlines():
                normalized = line.lower()
                if any(name in normalized for name in names):
                    values = re.findall(
                        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
                        line.replace(",", ""),
                    )
                    if values:
                        metrics[metric] = float(values[-1])
                        break

    if metrics["ipc"] == "" and metrics["instructions"] != "" and metrics["cycles"] not in ("", 0):
        metrics["ipc"] = metrics["instructions"] / metrics["cycles"]
    return metrics


def call_binary(
    cmd,
    timeout=None,
    perf_events=None,
    cpu_set=None,
    affinity_mask=None,
    vtune_result_dir=None,
    vtune_executable=None,
    vtune_bin=None,
    vtune_analysis="uarch-exploration",
):
    start = time.perf_counter()
    run_cmd = cmd
    if perf_events:
        run_cmd = ["perf", "stat", "-x,", "-e", PERF_EVENTS, "--", *cmd]
    if cpu_set:
        run_cmd = ["taskset", "-c", cpu_set, *run_cmd]
    if vtune_result_dir:
        run_cmd = [
            vtune_executable or "vtune", "-collect", vtune_analysis,
            "-result-dir", str(vtune_result_dir),
            "--", *run_cmd,
        ]
    child_env = os.environ.copy()
    if vtune_bin:
        child_env["PATH"] = str(vtune_bin) + os.pathsep + child_env.get("PATH", "")

    if affinity_mask is not None and os.name != "nt":
        raise RuntimeError("Windows affinity masks can only be applied on Windows")
    proc = subprocess.Popen(
        run_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=child_env,
    )
    if affinity_mask is not None:
        try:
            set_windows_process_affinity(proc, affinity_mask)
        except OSError:
            proc.kill()
            proc.communicate()
            raise

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
            if timeout is not None and time.perf_counter() - start > timeout:
                proc.kill()
                out, err = proc.communicate()
                raise TimeoutError(f"Command timed out after {timeout}s: {run_cmd}")
            if p:
                try:
                    processes = [p, *p.children(recursive=True)]
                    rss = sum(child.memory_info().rss for child in processes)
                    if peak_kb is None or rss > peak_kb:
                        peak_kb = rss
                except psutil.Error:
                    pass
            time.sleep(0.01)

        # last memory check
        if p:
            try:
                processes = [p, *p.children(recursive=True)]
                rss = sum(child.memory_info().rss for child in processes)
                if peak_kb is None or rss > peak_kb:
                    peak_kb = rss
            except psutil.Error:
                pass

        out, err = proc.communicate(timeout=1)

    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
        raise TimeoutError(f"Command timed out after {timeout}s: {run_cmd}")

    end = time.perf_counter()

    # convert memory
    if peak_kb is not None:
        peak_kb = peak_kb // 1024

    # decode
    out = out.decode("utf-8", errors="replace")
    err = err.decode("utf-8", errors="replace")

    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {proc.returncode}: {run_cmd}\n{err.strip()}"
        )
    return end - start, peak_kb, out.strip(), err.strip(), parse_perf_metrics(err)


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
        "binary","length","identity","threads","schedule_chunk","affinity","cpu_set",
        "affinity_mask","vtune_result_dir","vtune_report_path","repeat",
        "pairs_file","n_pairs","elapsed_s","peak_rss_kb",
        "l1_dcache_load_misses","instructions","cycles","ipc","stdout","stderr"
    ]

    raw_rows = []
    if args.perf and args.vtune:
        raise RuntimeError("Choose either --perf or --vtune; do not enable both")
    affinity_profiles = []
    if os.name == "nt":
        windows_profiles = [
            ("p_cores", args.p_core_mask),
            ("e_cores", args.e_core_mask),
            ("all_cores", args.all_core_mask),
        ]
        for name, mask_text in windows_profiles:
            if mask_text:
                try:
                    mask = int(mask_text, 0)
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid {name} affinity mask {mask_text!r}; use a hex mask such as 0x3"
                    ) from exc
                if mask <= 0:
                    raise ValueError(f"{name} affinity mask must be non-zero")
                affinity_profiles.append((name, "", mask))
    else:
        linux_profiles = [
            ("p_cores", args.p_core_cpus),
            ("e_cores", args.e_core_cpus),
            ("all_cores", args.all_core_cpus),
        ]
        for name, cpus in linux_profiles:
            if cpus:
                affinity_profiles.append((name, cpus, None))
    if affinity_profiles and os.name == "posix" and shutil.which("taskset") is None:
        raise RuntimeError("CPU affinity was requested but taskset was not found")
    if affinity_profiles and os.name not in ("posix", "nt"):
        raise RuntimeError("CPU affinity profiles are supported on Linux and Windows")
    if os.name == "nt" and any(
        (args.p_core_cpus, args.e_core_cpus, args.all_core_cpus)
    ):
        raise RuntimeError("Use --p-core-mask/--e-core-mask/--all-core-mask on Windows")
    if os.name == "posix" and any(
        (args.p_core_mask, args.e_core_mask, args.all_core_mask)
    ):
        raise RuntimeError("Use --p-core-cpus/--e-core-cpus/--all-core-cpus on Linux")
    if args.perf and (os.name != "posix" or shutil.which("perf") is None):
        raise RuntimeError("perf profiling was requested but Linux perf is unavailable")
    vtune_bin = Path(args.vtune_bin).expanduser().resolve() if args.vtune_bin else None
    if vtune_bin and not vtune_bin.is_dir():
        raise FileNotFoundError(f"VTune bin directory not found: {vtune_bin}")
    vtune_executable = None
    if vtune_bin:
        for candidate in ("vtune.exe", "vtune"):
            if (vtune_bin / candidate).is_file():
                vtune_executable = str(vtune_bin / candidate)
                break
        if not vtune_executable:
            raise FileNotFoundError(f"VTune CLI not found in {vtune_bin}")
    else:
        vtune_executable = shutil.which("vtune") or shutil.which("vtune.exe")
    if args.vtune and not vtune_executable:
        raise RuntimeError(
            "VTune was requested, but vtune was not found. Add its bin64 directory "
            "to PATH or pass --vtune-bin <path-to-bin64>."
        )
    if not affinity_profiles:
        affinity_profiles = [("unbound", "", None)]

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
                    for chunk_size in args.schedule_chunks:
                        for r in range(args.repeats):
                            for binary in args.binaries:
                                bin_path = Path(binary)
                                if not bin_path.is_absolute():
                                    bin_path = (Path(args.root_dir) / binary).resolve()

                                if not bin_path.exists():
                                    raise FileNotFoundError(f"Binary not found: {bin_path}")

                                if "tiled" in bin_path.name.lower():
                                    cmd = [str(bin_path), str(fasta_path), str(t), str(chunk_size)]
                                else:
                                    cmd = [
                                        str(bin_path), str(fasta_path),
                                        str(args.band_threshold), str(t), str(chunk_size)
                                    ]

                                for affinity, cpu_set, affinity_mask in affinity_profiles:
                                    vtune_result_dir = ""
                                    vtune_report_path = ""
                                    if args.vtune:
                                        result_name = (
                                            f"{bin_path.stem}_L{L}_id{int(ident * 100)}"
                                            f"_t{t}_chunk{chunk_size}_{affinity}_r{r}"
                                        )
                                        vtune_result_dir = str(
                                            outdir / "vtune_results" / result_name
                                        )
                                        vtune_report_path = str(
                                            Path(vtune_result_dir) / "summary.csv"
                                        )
                                    print(
                                        "[run]", affinity, cmd,
                                        f"cpu_set={cpu_set or '-'}",
                                        f"affinity_mask={hex(affinity_mask) if affinity_mask else '-'}",
                                        f"vtune_result={vtune_result_dir or '-'}",
                                    )

                                    elapsed, peak, out, err, metrics = call_binary(
                                        cmd,
                                        timeout=args.timeout,
                                        perf_events=PERF_EVENTS if args.perf else None,
                                        cpu_set=cpu_set or None,
                                        affinity_mask=affinity_mask,
                                        vtune_result_dir=vtune_result_dir or None,
                                        vtune_executable=vtune_executable,
                                        vtune_bin=vtune_bin,
                                        vtune_analysis=args.vtune_analysis,
                                    )
                                    if args.vtune:
                                        report_command = [
                                            vtune_executable,
                                            "-report", "summary",
                                            "-result-dir", vtune_result_dir,
                                            "-format", "csv",
                                            "-report-output", vtune_report_path,
                                        ]
                                        report_env = os.environ.copy()
                                        if vtune_bin:
                                            report_env["PATH"] = (
                                                str(vtune_bin) + os.pathsep
                                                + report_env.get("PATH", "")
                                            )
                                        try:
                                            subprocess.run(
                                                report_command,
                                                check=True,
                                                capture_output=True,
                                                text=True,
                                                timeout=args.timeout,
                                                env=report_env,
                                            )
                                        except subprocess.CalledProcessError as exc:
                                            raise RuntimeError(
                                                "VTune summary export failed: "
                                                f"{exc.stderr.strip()}"
                                            ) from exc
                                        vtune_metrics = parse_vtune_report(vtune_report_path)
                                        metrics.update({
                                            key: value for key, value in vtune_metrics.items()
                                            if value != ""
                                        })
                                    if args.perf and any(
                                        metrics[event] == ""
                                        for event in (
                                            "l1_dcache_load_misses", "instructions", "cycles"
                                        )
                                    ):
                                        print(
                                            "[warn] one or more perf counters are missing; "
                                            "check event support and perf permissions"
                                        )

                                    row = {
                                        "binary": bin_path.name,
                                        "length": L,
                                        "identity": ident,
                                        "threads": t,
                                        "schedule_chunk": chunk_size,
                                        "affinity": affinity,
                                        "cpu_set": cpu_set,
                                        "affinity_mask": (
                                            hex(affinity_mask) if affinity_mask is not None else ""
                                        ),
                                        "vtune_result_dir": vtune_result_dir,
                                        "vtune_report_path": vtune_report_path,
                                        "repeat": r,
                                        "pairs_file": str(fasta_path),
                                        "n_pairs": args.pairs_per_identity,
                                        "elapsed_s": elapsed,
                                        "peak_rss_kb": peak if peak else "",
                                        **metrics,
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
        key = (row["binary"], row["length"], row["identity"], row["threads"],
               row["schedule_chunk"], row["affinity"], row["cpu_set"],
               row["affinity_mask"])
        groups.setdefault(key, []).append(row)

    with agg_csv.open("w", newline="", encoding="utf-8") as af:
        fn2 = ["binary","length","identity","threads","schedule_chunk","affinity","cpu_set",
               "affinity_mask",
               "n_runs","median_time","mean_time","stdev_time","median_rss_kb",
               "median_l1_dcache_load_misses","median_instructions","median_cycles","median_ipc"]
        w2 = csv.DictWriter(af, fieldnames=fn2)
        w2.writeheader()

        for key, rows in sorted(groups.items()):
            times = [r["elapsed_s"] for r in rows]
            rss = [int(r["peak_rss_kb"]) for r in rows if r["peak_rss_kb"]]
            counter_columns = (
                "l1_dcache_load_misses", "instructions", "cycles", "ipc"
            )
            rec = {
                "binary": key[0],
                "length": key[1],
                "identity": key[2],
                "threads": key[3],
                "schedule_chunk": key[4],
                "affinity": key[5],
                "cpu_set": key[6],
                "affinity_mask": key[7],
                "n_runs": len(rows),
                "median_time": median(times),
                "mean_time": mean(times),
                "stdev_time": stdev(times) if len(times) > 1 else 0.0,
                "median_rss_kb": median(rss) if rss else "",
            }
            for column in counter_columns:
                values = [r[column] for r in rows if r[column] != ""]
                rec[f"median_{column}"] = median(values) if values else ""
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
        for (bname, affinity, chunk), g in tmp.groupby(
            ["binary", "affinity", "schedule_chunk"]
        ):
            g2 = g.sort_values("threads")
            label = f"{bname}, {affinity}, chunk={chunk}"
            plt.plot(g2.threads, g2.median_time, marker='o', label=label)
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

        speedup_fig = plt.figure(figsize=(7,4))
        for (bname, affinity, chunk), g in tmp.groupby(
            ["binary", "affinity", "schedule_chunk"]
        ):
            baseline = g.loc[g.threads == 1, "median_time"]
            if baseline.empty or baseline.iloc[0] <= 0:
                continue
            g2 = g.sort_values("threads")
            speedups = baseline.iloc[0] / g2.median_time
            plt.plot(
                g2.threads, speedups, marker='o',
                label=f"{bname}, {affinity}, chunk={chunk}"
            )
        max_threads = int(tmp.threads.max())
        plt.plot([1, max_threads], [1, max_threads], linestyle="--", color="gray",
                 label="Ideal linear speedup")
        plt.xlabel("Threads")
        plt.ylabel("Speedup vs. one thread")
        plt.title(f"Speedup at length={Lmax}, identity={int(ident*100)}%")
        plt.grid(True)
        plt.legend()
        speedup_fig.tight_layout()
        speedup_path = outdir / f"speedup_L{Lmax}_id{int(ident*100)}.png"
        speedup_fig.savefig(speedup_path, dpi=200)
        plt.close(speedup_fig)
        print("[done] Plot:", speedup_path)

# ============================================================
# Auto-compile (optional)
# ============================================================
def auto_compile(args):
    srcdir = Path(args.root_dir) / "src"
    tiled = srcdir / "tiled_sw.c"
    wfa = srcdir / "wfa_banded.c"
    compiler = shutil.which(args.compiler) if args.compiler else None
    if args.compiler and compiler is None:
        raise FileNotFoundError(f"C compiler not found: {args.compiler}")
    compiler_candidates = ("gcc", "cl", "clang") if os.name == "nt" else ("gcc", "clang", "cl")
    if compiler is None:
        compiler = next(
            (shutil.which(candidate) for candidate in compiler_candidates
             if shutil.which(candidate)),
            None,
        )
    if compiler is None:
        raise RuntimeError("No supported C compiler found (tried gcc, clang, and cl)")

    compiler_name = Path(compiler).name.lower()
    msvc = compiler_name in ("cl", "cl.exe")
    sources = (
        (tiled, Path(args.root_dir) / "tiled_sw.exe"),
        (wfa, Path(args.root_dir) / "wfa_banded.exe"),
    )
    for source, output in sources:
        if not source.exists():
            raise FileNotFoundError(f"Source file not found: {source}")
        print("[compile]", source.name, "with", compiler_name)
        if msvc:
            with tempfile.TemporaryDirectory(prefix="alignment-build-") as obj_dir:
                subprocess.run(
                    [
                        compiler, "/O2", "/openmp", f"/Fo{obj_dir}{os.sep}",
                        f"/Fe{output}", str(source),
                    ],
                    check=True,
                    cwd=args.root_dir,
                )
        else:
            subprocess.run(
                [compiler, "-O3", "-fopenmp", str(source), "-o", str(output)],
                check=True,
                cwd=args.root_dir,
            )



# ============================================================
# CLI
# ============================================================
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default=str(PROJECT_ROOT / "results"))
    p.add_argument("--pairs-per-identity", type=int, default=200)
    p.add_argument("--lengths", nargs="+", type=int, default=[2000,5000])
    p.add_argument("--identities", nargs="+", type=float, default=[0.98,0.95,0.90,0.80])
    p.add_argument("--threads-list", nargs="+", type=int, default=[1,2,4,8])
    p.add_argument("--schedule-chunks", nargs="+", type=int, default=[1])
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--band-threshold", type=int, default=200)
    p.add_argument("--perf", action="store_true",
                   help="Wrap runs with Linux perf stat")
    p.add_argument("--vtune", action="store_true",
                   help="Collect VTune uarch-exploration data and export a CSV summary")
    p.add_argument("--vtune-analysis", default="uarch-exploration",
                   help="VTune collection type (default: uarch-exploration)")
    p.add_argument("--vtune-bin",
                   help="Path to VTune's bin64 directory; added to child-process PATH")
    p.add_argument("--p-core-cpus", help="Linux taskset CPU list for P-cores, e.g. 0-5")
    p.add_argument("--e-core-cpus", help="Linux taskset CPU list for E-cores, e.g. 6-11")
    p.add_argument("--all-core-cpus", help="Linux taskset CPU list for all selected cores")
    p.add_argument("--p-core-mask", help="Windows CPU affinity mask for P-cores, e.g. 0x3")
    p.add_argument("--e-core-mask", help="Windows CPU affinity mask for E-cores, e.g. 0xC")
    p.add_argument("--all-core-mask", help="Windows CPU affinity mask for all selected cores")
    p.add_argument("--binaries", nargs="+",
                  default=["tiled_sw.exe", "wfa_banded.exe"])
    p.add_argument("--root-dir", default=str(PROJECT_ROOT))
    p.add_argument("--compiler", help="C compiler for --auto-compile (gcc, clang, or cl)")
    p.add_argument("--auto-compile", action="store_true")

    args = p.parse_args()
    if any(value < 1 for value in args.threads_list + args.schedule_chunks):
        p.error("thread counts and schedule chunks must be positive")

    if args.auto_compile:
        auto_compile(args)

    run_sweep(args)


if __name__ == "__main__":
    main()
