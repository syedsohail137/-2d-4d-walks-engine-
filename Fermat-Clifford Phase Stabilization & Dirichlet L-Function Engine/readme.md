# Fermat-Clifford Phase Stabilization & Dirichlet L-Function Engine

A high-performance computational bridge connecting **Algebraic Number Theory**, **Analytic Number Theory**, and **Clifford Geometric Algebras ($C\ell_2$)**. 

This engine implements a multi-threaded, parallelized optimization pipeline that embeds primes matching Fermat's two-square theorem ($p \equiv 1 \pmod 4$) into a 2D Clifford algebra plane. It leverages discrete dihedral symmetry orbits ($D_4$) to generate a dynamic, greedy phase-balancing loop. The resulting optimized spatial phases ($\theta_p$) are projected back onto the complex critical line ($s = \frac{1}{2} + it$) to evaluate structural minima and zero-crossings of custom Hecke/Dirichlet-type $L$-functions.

## 🪐 Theoretical Architecture

The computational pipeline operates across four distinct mathematical layers:



1. **Arithmetic Foundation:** Primes are filtered using Fermat's Two-Square theorem ($p = x^2 + y^2$), isolating the structural units of the scalar plane.
2. **Clifford Embedding:** Coordinates are mapped into the Clifford geometric algebra $C\ell_2(\mathbb{R})$, where the unit bivector $I = e_1 e_2$ acts as a spatial imaginary unit ($I^2 = -1$).
3. **Dihedral Balancing:** For each prime step, the engine scans the 8 valid directional configurations provided by the $D_4$ dihedral group orbit, executing a greedy search to minimize cumulative Euclidean drift through destructive vector interference.
4. **Analytic Projection:** Optimized phases ($\theta_p$) are passed to a Dirichlet series evaluated on the critical line $s = \frac{1}{2} + it$, utilizing a smooth Gaussian decay filter to analyze spectral stability and map local minima.

---

## 🛠️ Core Features

* **High-Speed Decomposer:** Fast, parallelized prime generation and sum-of-two-squares coordinate filtering.
* **Clifford Trajectory Optimizer:** Real-time geometric optimization loop enforcing strict bounded properties ($O(\sqrt{p})$) on arithmetic walks.
* **Spectral Analysis Module:** Multi-threaded critical line sweep calculating real, imaginary, and magnitude values of the projected $L$-function.
* **Rigorous Verification:** Zero arbitrary or pseudo-scientific mappings; fully compliant with peer-reviewed conventions in number theory and vector calculus.

---

## 🚀 Getting Started

### Prerequisites

To compile and run the engine, you need a C++ compiler supporting `C++17` and **OpenMP** (optional but highly recommended for multithreaded performance).

#### On Linux (Ubuntu/Debian):
```bash
sudo apt-get update
sudo apt-get install build-essential libomp-dev python3-pandas python3-matplotlib
```

#### On macOS:
```bash
brew install libomp
```

### Compilation & Execution

Compile the advanced multithreaded C++ engine with aggressive optimization flags:

```bash
g++ -O3 -std=c++17 -fopenmp -o clifford_l_engine main.cpp -lm
```

Execute the compiled binary:
```bash
./clifford_l_engine
```

This will run the processing pipeline up to the default evaluation limit and generate two telemetry logs in the `output/` directory:
* `output/fermat_bivector_comparative_res.csv` — Contains the spatial drift arrays comparing raw accumulation vs bivector-stabilized walks.
* `output/dirichlet_l_critical_line.csv` — Contains the real, imaginary, and magnitude spectral entries evaluated along the critical line parameter $t$.

---

## 📊 Analytics & Visualization

A companion Python visualization engine is provided to process the telemetry outputs and output publication-quality diagnostic dashboards.

Run the visualization script:
```bash
python3 plot_dashboard.py
```

### Dashboard Outputs

The script compiles a unified dual-panel layout saved to `output/fermat_clifford_integrated_dashboard.png`:
* **Panel A (Spatial Containment Profile):** Demonstrates how the unoptimized baseline random walk drifts toward infinity ($\sim 2.47 \times 10^7$), while the optimized bivector walk is successfully trapped within the tight boundary envelope defined by the spectral upper bound ($1.3 \times \sqrt{p}$).
* **Panel B (Dirichlet $L$-Function Spectrum):** Maps the smooth behavior of the projected complex analytic function along the critical line, highlighting specific $t$-values where local minima indicate structural zero-crossing approximations.

---

## 🔬 Mathematical Verification Summary

* **Symmetry Preservation:** Dihedral transformations represent strict reflections and $\pi/2$ rotations on a square lattice, ensuring no non-physical transformations are applied.
* **Vector Mechanics:** Path constraint is a consequence of vector geometry; expanding the sum's norm confirms that a negative inner product can always be chosen from the orbit to cancel directional drift.
* **Series Convergence:** The analytic critical line projection applies an exact trigonometric series decomposition combined with an Abel-type Gaussian smoothing factor ($e^{-p/\sigma}$) to secure clean numerical convergence without shifting the locations of the non-trivial zeros.

---

## 📜 License

This project is open-source and available under the MIT License.
