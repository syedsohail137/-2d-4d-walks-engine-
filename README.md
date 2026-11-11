# 2D–4D Walks Engine

This repository collects several standalone C++ experiments related to bounded random walks, prime lattice mappings, stochastic simulation, and spectral/number-theoretic analysis. The modules are intentionally independent and are organized as separate subdirectories.

## Repository layout

- `src/` — root demo; builds a small C++17 program that generates primes, maps them into four-square decompositions, and writes CSV output.
- `v6reverse/` — reverse lookup engine for the first Riemann zeros and prime-based determinant search.
- `b4-engine-4d-walks-main/` — 4D walk engine project with multiple variants (`v2`, `V3`, etc.).
- `Time-Variant Bounded Random Walk/` — bounded walk and discrepancy-minimization experiments.
- `Bound-GBM-CIR/` and `Bound-GBM-CIR-MonteCarlo/` — SDE simulation and lattice-bound validation studies.
- `Discrete Lyapunov Optimization/` — Lyapunov-style control experiments over hyperoctahedral transformations.
- `Fermat-Clifford Phase Stabilization & Dirichlet L-Function Engine/` — prime filtering with Clifford and Dirichlet-style spectral analysis.
- `SO2SO3Clifford/` — Clifford and rotation-group experimental work.
- `asymptotic-hyperspace/` — quaternionic / hyperspace prototype and supporting documentation.
- `comparison/` — comparison utilities and output files.

## Root demo build

Requirements:
- CMake 3.10+
- C++17 compiler

Build:

```bash
cmake -S src -B build
cmake --build build --config Release
```

Run:

```bash
./build/4D_Engine
```

The program writes CSV output in the current working directory and prints summary information to the terminal.

## Notes

- Each subdirectory is a separate project or research artifact.
- Some directories contain experimental scripts, PDFs, and generated output files rather than a conventional build.
- For precise usage for a given subproject, consult the README in that directory.

## License

MIT License.
