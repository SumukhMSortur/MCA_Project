import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================
# Load aggregated CSV
# ==========================
df = pd.read_csv("aggregated_results.csv")

sns.set(style="whitegrid")


# ==========================================================
# 1) MEMORY (RSS) VS THREADS — OVERALL FOR EACH LENGTH
# ==========================================================
def plot_rss_vs_threads_overall():
    lengths = sorted(df["length"].unique())

    for L in lengths:
        subset = df[df["length"] == L]

        plt.figure(figsize=(12, 7))
        sns.lineplot(
            data=subset,
            x="threads",
            y="median_rss_kb",
            hue="binary",
            marker="o"
        )

        plt.title(f"Memory (RSS) vs Threads (length = {L})")
        plt.xlabel("Threads")
        plt.ylabel("Median RSS (KB)")
        plt.grid(True, alpha=0.3)

        out = f"rss_vs_threads_length{L}.png"
        plt.savefig(out, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[✔] Saved: {out}")


# ==========================================================
# 2) MEMORY (RSS) VS THREADS — PER IDENTITY PER LENGTH
# ==========================================================
def plot_rss_vs_threads_for_all_identities():
    lengths = sorted(df["length"].unique())
    identities = sorted(df["identity"].unique())

    for L in lengths:
        for ID in identities:
            subset = df[(df["length"] == L) & (df["identity"] == ID)]
            if subset.empty:
                continue

            plt.figure(figsize=(12, 7))
            sns.lineplot(
                data=subset,
                x="threads",
                y="median_rss_kb",
                hue="binary",
                marker="o"
            )

            plt.title(f"Memory (RSS) vs Threads (length = {L}, identity = {ID:.2f})")
            plt.xlabel("Threads")
            plt.ylabel("Median RSS (KB)")
            plt.grid(True, alpha=0.3)

            out = f"rss_vs_threads_length{L}_id{int(ID * 100)}.png"
            plt.savefig(out, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"[✔] Saved: {out}")


# ==========================================================
# 3) RUNTIME VS IDENTITY — OVERALL FOR EACH THREAD COUNT
# ==========================================================
def plot_runtime_vs_identity_overall():
    threads_list = sorted(df["threads"].unique())

    for T in threads_list:
        subset = df[df["threads"] == T]

        plt.figure(figsize=(12, 7))
        sns.lineplot(
            data=subset,
            x="identity",
            y="median_time",
            hue="binary",
            marker="o"
        )

        plt.title(f"Runtime vs Identity (threads = {T})")
        plt.xlabel("Identity (%)")
        plt.ylabel("Median Runtime (s)")
        plt.grid(True, alpha=0.3)

        out = f"runtime_vs_identity_threads{T}.png"
        plt.savefig(out, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[✔] Saved: {out}")


# ==========================================================
# 4) RUNTIME VS IDENTITY — PER LENGTH SEPARATELY (OPTIONAL)
# ==========================================================
def plot_runtime_vs_identity_per_length():
    lengths = sorted(df["length"].unique())

    for L in lengths:
        subset = df[df["length"] == L]

        plt.figure(figsize=(12, 7))
        sns.lineplot(
            data=subset,
            x="identity",
            y="median_time",
            hue="binary",
            marker="o"
        )

        plt.title(f"Runtime vs Identity (length = {L})")
        plt.xlabel("Identity (%)")
        plt.ylabel("Median Runtime (s)")
        plt.grid(True, alpha=0.3)

        out = f"runtime_vs_identity_length{L}.png"
        plt.savefig(out, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[✔] Saved: {out}")


# ==========================================================
# 5) SPEEDUP BY AFFINITY PROFILE AND DYNAMIC CHUNK
# ==========================================================
def plot_affinity_speedup():
    required = {"affinity", "schedule_chunk"}
    if not required.issubset(df.columns):
        print("[skip] Affinity speedup columns are not present; rerun the updated harness")
        return

    lengths = sorted(df["length"].unique())
    identities = sorted(df["identity"].unique())

    for length in lengths:
        for identity in identities:
            subset = df[(df["length"] == length) & (df["identity"] == identity)]
            if subset.empty:
                continue
            plt.figure(figsize=(12, 7))

            for (binary, affinity, chunk), group in subset.groupby(
                ["binary", "affinity", "schedule_chunk"]
            ):
                baseline = group.loc[group["threads"] == 1, "median_time"]
                if baseline.empty or baseline.iloc[0] <= 0:
                    continue
                group = group.sort_values("threads")
                speedup = baseline.iloc[0] / group["median_time"]
                plt.plot(
                    group["threads"], speedup, marker="o",
                    label=f"{binary}, {affinity}, chunk={chunk}"
                )

            max_threads = int(subset["threads"].max())
            plt.plot(
                [1, max_threads], [1, max_threads], linestyle="--",
                color="gray", label="Ideal linear speedup"
            )
            plt.title(
                f"Batch alignment speedup (length={length}, identity={identity:.2f})"
            )
            plt.xlabel("Threads")
            plt.ylabel("Speedup vs. one thread")
            plt.grid(True, alpha=0.3)
            plt.legend()

            out = f"speedup_affinity_length{length}_id{int(identity * 100)}.png"
            plt.savefig(out, dpi=300, bbox_inches="tight")
            plt.close()
            print(f"[✔] Saved: {out}")


# ==========================================================
# 6) PERF METRICS — IPC AND L1 DATA-CACHE MISSES
# ==========================================================
def plot_perf_metrics():
    metric_columns = {"median_ipc", "median_l1_dcache_load_misses"}
    if not metric_columns.issubset(df.columns):
        print("[skip] Hardware counters are not present; run the harness with --perf")
        return

    if df["median_ipc"].notna().any():
        plt.figure(figsize=(12, 7))
        sns.lineplot(
            data=df,
            x="threads",
            y="median_ipc",
            hue="affinity" if "affinity" in df.columns else "binary",
            style="binary",
            marker="o"
        )
        plt.title("Instruction per cycle by thread count")
        plt.xlabel("Threads")
        plt.ylabel("Median IPC")
        out = "ipc_vs_threads.png"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"[✔] Saved: {out}")

    if df["median_l1_dcache_load_misses"].notna().any():
        plt.figure(figsize=(12, 7))
        sns.lineplot(
            data=df,
            x="threads",
            y="median_l1_dcache_load_misses",
            hue="affinity" if "affinity" in df.columns else "binary",
            style="binary",
            marker="o"
        )
        plt.title("L1 data-cache load misses by thread count")
        plt.xlabel("Threads")
        plt.ylabel("Median L1 data-cache load misses")
        out = "l1_dcache_misses_vs_threads.png"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"[✔] Saved: {out}")


# ==========================================================
# MAIN FUNCTION — RUN EVERYTHING
# ==========================================================
if __name__ == "__main__":
    print("\n=== Generating All Graphs ===\n")

    plot_rss_vs_threads_overall()
    plot_rss_vs_threads_for_all_identities()

    plot_runtime_vs_identity_overall()
    plot_runtime_vs_identity_per_length()
    plot_affinity_speedup()
    plot_perf_metrics()

    print("\n=== All plots generated successfully! ===\n")
    plot_runtime_vs_identity_overall()
    plot_runtime_vs_identity_per_length()

    print("\n=== All plots generated successfully! ===\n")










