import os
import pandas as pd
import matplotlib.pyplot as plt

def render_stabilisation_dashboard(walk_csv="output/fermat_bivector_comparative_res.csv", 
                                   spectral_csv="output/dirichlet_l_critical_line.csv"):
    """
    Generates a high-resolution dual-panel diagnostic dashboard analyzing
    the spatial field containment profile and the spectral critical line behavior.
    """
    if not os.path.exists(walk_csv) or not os.path.exists(spectral_csv):
        print("[ERROR] Analytics matrices missing from the output directory.")
        print("[ERROR] Please run your multithreaded C++ engine first.")
        return

    print("[SYSTEM] Importing telemetric logs into pandas parser...")
    df_walk = pd.read_csv(walk_csv)
    df_spec = pd.read_csv(spectral_csv)

    # Initialize master layout canvas
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 7), dpi=300)
    plt.subplots_adjust(wspace=0.25)

    # --- PANEL 1: SPATIAL FIELD CONTAINMENT ---
    color_a = '#d95f02'
    color_b = '#1f77b4'
    color_bound = '#2ca02c'

    ax1_twin = ax1.twinx()
    
    # Plot Track A (Raw Drift) on the left Y-axis
    ax1.set_xlabel('Fermat Prime Sequence (p)', fontsize=11, fontweight='bold', labelpad=8)
    ax1.set_ylabel('Track A: Unoptimized Drift Scale', color=color_a, fontsize=11, fontweight='bold')
    line1 = ax1.plot(df_walk['Prime'], df_walk['Track_A_Raw'], color=color_a, linestyle='--', alpha=0.7, label='Track A: Raw Unoptimized Drift')
    ax1.tick_params(axis='y', labelcolor=color_a)
    ax1.grid(True, linestyle=':', alpha=0.5)

    # Plot Track B (Greedy Balanced Walk) and Boundary on the right Y-axis
    ax1_twin.set_ylabel('Track B & Bound Scale Magnitude', color=color_b, fontsize=11, fontweight='bold')
    line2 = ax1_twin.plot(df_walk['Prime'], df_walk['Track_B_Greedy'], color=color_b, linewidth=1.5, label='Track B: Bivector Balanced Walk')
    line3 = ax1_twin.plot(df_walk['Prime'], df_walk['Spectral_Bound'], color=color_bound, linewidth=1.75, linestyle=':', label='Spectral Boundary (1.3 * sqrt(p))')
    ax1_twin.tick_params(axis='y', labelcolor=color_b)

    # Conjoin legends
    lines_1 = line1 + line2 + line3
    labels_1 = [ln.get_label() for ln in lines_1]
    ax1.legend(lines_1, labels_1, loc='upper left', frameon=True, facecolor='white', edgecolor='#e0e0e0')
    ax1.set_title('A: Coordinate Drift Space Containment Profile', fontsize=12, fontweight='bold', pad=10)

    # --- PANEL 2: SPECTRUM & ZERO CROSSINGS ---
    # Find local minima locations to flag potential zero crossings
    minima_df = df_spec[df_spec['Magnitude'] == df_spec['Magnitude'].min()]
    
    ax2.plot(df_spec['t'], df_spec['Magnitude'], color='#7570b3', linewidth=1.75, label=r'|L($\frac{1}{2}$ + it)|')
    ax2.plot(df_spec['t'], df_spec['Real'], color='#e7298a', linewidth=1.0, linestyle='-.', alpha=0.6, label='Real Component')
    ax2.plot(df_spec['t'], df_spec['Imag'], color='#66a61e', linewidth=1.0, linestyle='-.', alpha=0.6, label='Imaginary Component')
    
    # Highlight the absolute minimum value observed in the telemetry logs
    for idx, row in minima_df.iterrows():
        ax2.scatter(row['t'], row['Magnitude'], color='red', s=40, zorder=5)
        ax2.annotate(f"Min |L|={row['Magnitude']:.4f}\nat t={row['t']:.2f}",
                     xy=(row['t'], row['Magnitude']),
                     xytext=(row['t'] + 5, row['Magnitude'] + 0.3),
                     arrowprops=dict(facecolor='black', arrowstyle='->', linestyle=':'),
                     fontsize=9, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.3, ec='orange'))

    ax2.set_xlabel('Spectral Axis Parameter (t)', fontsize=11, fontweight='bold', labelpad=8)
    ax2.set_ylabel('Analytic Functional Value Scale', fontsize=11, fontweight='bold', labelpad=8)
    ax2.axhline(0, color='black', linewidth=0.8, linestyle='-')
    ax2.grid(True, linestyle=':', alpha=0.5)
    ax2.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='#e0e0e0')
    ax2.set_title(r'B: Dirichlet L-Function Behavior along Critical Line ($s=\frac{1}{2}+it$)', fontsize=12, fontweight='bold', pad=10)

    # Output master rendering
    master_title = "Fermat-Clifford Integrated Analytics Dashboard\nSymmetries of the Dual-Basis Scalar Plane"
    plt.suptitle(master_title, fontsize=15, fontweight='bold', y=1.02)
    
    output_png = "output/fermat_clifford_integrated_dashboard.png"
    plt.savefig(output_png, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"[COMPLETE] Dashboard graphic successfully output to: {output_png}")

if __name__ == "__main__":
    render_stabilisation_dashboard()
