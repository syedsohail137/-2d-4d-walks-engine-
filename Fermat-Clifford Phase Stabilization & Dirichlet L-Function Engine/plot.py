import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_fermat_clifford_telemetry(csv_path="output/fermat_bivector_comparative_res.csv"):
    """
    Parses telemetric logging output files and plots Track A (Raw Drift) 
    and Track B (Bivector Balanced Walk) against the Spectral Upper Bound.
    """
    if not os.path.exists(csv_path):
        print(f"[ERROR] Telemetry file not detected at path: {csv_path}")
        print("[ERROR] Please execute the optimized C++ binary compilation loop first.")
        return

    print(f"[SYSTEM] Parsing telemetric data frame from: {csv_path}")
    
    # Efficient parsing chunk allocation via pandas engine
    df = pd.read_csv(csv_path)
    
    # Initialize structural plot canvas
    plt.figure(figsize=(12, 7), dpi=300)
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # Plot Track A (Raw Drift) on a dual axis if scales differ immensely
    fig, ax1 = plt.subplots(figsize=(12, 7))
    
    color_track_a = '#e41a1c'
    color_track_b = '#377eb8'
    color_bound   = '#4daf4a'
    
    ax1.set_xlabel('Fermat Prime Value (p)', fontsize=12, fontweight='bold', labelpad=10)
    ax1.set_ylabel('Track A: Unoptimized Raw Drift Magnitude', color=color_track_a, fontsize=12, fontweight='bold')
    line1 = ax1.plot(df['Prime'], df['Track_A_Raw'], color=color_track_a, linestyle='--', alpha=0.8, label='Track A: Raw Unoptimized Drift')
    ax1.tick_params(axis='y', labelcolor=color_track_a)
    ax1.grid(True, which='both', linestyle=':', alpha=0.5)
    
    # Instantiate secondary shared-X axis for high-resolution tracking of Track B and Bounds
    ax2 = ax1.twinx()  
    ax2.set_ylabel('Track B & Spectral Bound Scale Magnitude', color=color_track_b, fontsize=12, fontweight='bold')
    
    line2 = ax2.plot(df['Prime'], df['Track_B_Greedy'], color=color_track_b, linewidth=1.75, label='Track B: Bivector Balanced Walk')
    line3 = ax2.plot(df['Prime'], df['Spectral_Bound'], color=color_bound, linewidth=2.0, linestyle=':', label='Target Boundary Envelope (1.3 * sqrt(p))')
    ax2.tick_params(axis='y', labelcolor=color_track_b)
    
    # Consolidate legend arrays across dual axes matrices
    lines = line1 + line2 + line3
    labels = [ln.get_label() for ln in lines]
    ax1.legend(lines, labels, loc='upper left', frameon=True, facecolor='white', edgecolor='#d3d3d3', fontsize=10)
    
    plt.title('Asymptotic Field Containment Analysis:\nFermat-Clifford Stabilisation Loop Trajectories', fontsize=14, fontweight='bold', pad=15)
    
    # Save diagnostic visualization asset
    output_plot_path = "output/fermat_clifford_containment_profile.png"
    plt.savefig(output_plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"[COMPLETE] Diagnostic visualization plot compiled successfully to: {output_plot_path}")

if __name__ == "__main__":
    plot_fermat_clifford_telemetry()
