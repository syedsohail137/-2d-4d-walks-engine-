# 2D–4D Walks Engine

A collection of C++ experiments exploring bounded walks, prime-lattice coordinates, discrepancy minimization, and phase-modulated trajectories. The repository also contains supporting mathematical write-ups, generated data, plots, and Python visualization code.

> **Scope note:** This repository contains several related but separately organized experiments. The root-level build instructions below apply specifically to the small C++17 executable in `src/`; other folders may have their own code and instructions. The mathematical models and empirical results are exploratory and should not be read as formal proofs or production guarantees.

## Repository contents

| Path | Contents |
| --- | --- |
| `src/` | Root C++17 demo: generates primes up to 2,000, decomposes them into four-square coordinates, applies phase modulation and sign selection, and writes a CSV results file. |
| `b4-engine-4d-walks-main/` | A separate 4D walk simulation project, with its own README and source tree. |
| `Time-Variant Bounded Random Walk/` | 2D/4D bounded-walk experiments and supporting paper/data files. |
| `Bound-GBM-CIR/`, `Bound-GBM-CIR-MonteCarlo/` | Bounded stochastic-process experiments. |
| `Discrete Lyapunov Optimization/`, `Fermat-Clifford Phase Stabilization & Dirichlet L-Function Engine/`, `SO2SO3Clifford/`, `asymptotic-hyperspace/`, `comparison/` | Additional optimization, phase, geometry, and comparison experiments. |
| `paper_draft.tex` and PDF files | Mathematical notes and project write-ups. |
| `plot_generator.py`, `drift_compression_plot.png` | Plotting script and example plot. |

## Build and run the root demo

### Requirements

- CMake 3.10 or newer
- A C++17-capable compiler (for example, GCC, Clang, or MSVC)

### Build

Run these commands from the repository root:

```bash
cmake -S src -B build
cmake --build build --config Release
```

The executable is named `4D_Engine` (on Windows, it may be under `build/Release/`). For example, on Linux or macOS:

```bash
./build/4D_Engine
```

The program writes `mertens_drift_high_res.csv` to its **current working directory** and prints selected results to the terminal. The current root demo uses a prime limit of 2,000, embedded phase constants, and a square-root displacement comparison; these settings are defined in `src/main.cpp`.

## Working with the other experiments

The other directories are standalone projects or research artifacts rather than targets of the root `src/CMakeLists.txt`. See any directory-specific README for its own prerequisites and commands. Results may depend on the particular implementation and parameters used in that subproject.

## Research materials

The repository includes papers and LaTeX drafts discussing the mathematical motivation behind the experiments. These documents describe proposed models and analyses; validate assumptions and reproduce results independently before relying on them.

## License

No license is specified at the repository root. Unless a license is added, reuse and redistribution are subject to the repository owner's rights. Check individual subdirectories for any separate license files.
