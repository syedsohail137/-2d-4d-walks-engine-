import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else "spectral_comparison_large.csv"
    df = pd.read_csv(csv_path)

    n = len(df)
    idx = df["Index"].values
    primes = df["Prime"].values

    r_exact = (df["Riemann_State"] == "Exact Match").astype(int).values
    b_exact = (df["Base_State"] == "Exact Match").astype(int).values

    # Cumulative exact match percentages
    r_cum_pct = 100.0 * np.cumsum(r_exact) / idx
    b_cum_pct = 100.0 * np.cumsum(b_exact) / idx

    # Rolling Z-score of Riemann vs. Paired Baseline
    cum_diff = np.cumsum(r_exact - b_exact)
    pooled_p = np.clip((np.cumsum(r_exact) + np.cumsum(b_exact)) / (2.0 * idx), 1e-6, 1.0 - 1e-6)
    se = np.sqrt(2.0 * idx * pooled_p * (1.0 - pooled_p))
    rolling_z = cum_diff / se

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=150)
    fig.suptitle(
        f"Riemann-Prime Spectral Lattice Projection Diagnostics (N = {n:,}, p_max = {primes[-1]:,})",
        fontsize=14,
        fontweight="bold",
    )

    # Panel 1: Cumulative Exact Match Rate Decay
    ax1 = axes[0, 0]
    ax1.plot(idx, r_cum_pct, label="Riemann Zero Phase", color="#1f77b4", linewidth=2)
    ax1.plot(idx, b_cum_pct, label="Uniform Random Phase", color="#ff7f0e", linestyle="--", linewidth=1.8)
    ax1.set_xscale("log" if n > 200 else "linear")
    ax1.set_title("Cumulative Exact Match Rate (%) vs. Node Index")
    ax1.set_xlabel("Node Index n")
    ax1.set_ylabel("Cumulative Exact Match (%)")
    ax1.legend()

    # Panel 2: Rolling Z-Score Trajectory
    ax2 = axes[0, 1]
    ax2.plot(idx, rolling_z, color="#2ca02c", linewidth=1.8, label="Cumulative Z-Score (Riemann vs. Base)")
    ax2.axhline(0.0, color="black", linestyle="-", linewidth=1)
    ax2.axhline(1.96, color="#d62728", linestyle=":", linewidth=1.2, label="±1.96 (p = 0.05 threshold)")
    ax2.axhline(-1.96, color="#d62728", linestyle=":", linewidth=1.2)
    ax2.set_xscale("log" if n > 200 else "linear")
    ax2.set_title("Cumulative Exact-Match Z-Score Trajectory")
    ax2.set_xlabel("Node Index n")
    ax2.set_ylabel("Z-Score")
    ax2.legend()

    # Panel 3: Normalized Initial Determinant Distribution
    ax3 = axes[1, 0]
    norm_r_d0 = df["Riemann_D_Init"] / (df["Prime"] / 2.0)
    norm_b_d0 = df["Base_D_Init"] / (df["Prime"] / 2.0)
    bins = np.linspace(-1.0, 1.0, 41)
    ax3.hist(norm_r_d0, bins=bins, alpha=0.6, label="Riemann D_0 / (p/2)", color="#1f77b4", density=True)
    ax3.hist(norm_b_d0, bins=bins, alpha=0.5, label="Baseline D_0 / (p/2)", color="#ff7f0e", density=True)
    ax3.set_title("Normalized Initial Determinant Geometry")
    ax3.set_xlabel("Normalized Determinant D_0 / (p / 2)")
    ax3.set_ylabel("Probability Density")
    ax3.legend()

    # Panel 4: Search Steps Distribution Comparison
    ax4 = axes[1, 1]
    max_plot_step = int(min(max(df["Riemann_Steps"].max(), df["Base_Steps"].max()), 25))
    step_bins = np.arange(0, max_plot_step + 2) - 0.5
    ax4.hist(df["Riemann_Steps"], bins=step_bins, alpha=0.65, label="Riemann Search Steps", color="#1f77b4", density=True)
    ax4.hist(df["Base_Steps"], bins=step_bins, alpha=0.55, label="Baseline Search Steps", color="#ff7f0e", density=True)
    ax4.set_title("Lattice Convergence Search Radius Distribution")
    ax4.set_xlabel("Search Steps |D_conv - D_0|")
    ax4.set_ylabel("Fraction of Primes")
    ax4.legend()

    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    out_img = csv_path.rsplit(".", 1)[0] + "_diagnostics.png"
    plt.savefig(out_img)
    print(f"Saved diagnostic plot to {out_img}")


if __name__ == "__main__":
    main()
