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
# MAIN FUNCTION — RUN EVERYTHING
# ==========================================================
if __name__ == "__main__":
    print("\n=== Generating All Graphs ===\n")

    plot_rss_vs_threads_overall()
    plot_rss_vs_threads_for_all_identities()

    plot_runtime_vs_identity_overall()
    plot_runtime_vs_identity_per_length()

    print("\n=== All plots generated successfully! ===\n")
    plot_runtime_vs_identity_overall()
    plot_runtime_vs_identity_per_length()

    print("\n=== All plots generated successfully! ===\n")












