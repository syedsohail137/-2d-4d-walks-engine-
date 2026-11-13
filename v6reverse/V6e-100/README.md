# V6e-100: Coupled spectral lattice and phase-accumulator study

This directory contains a C++ prototype for testing a coupled spectral-lattice model against the first non-trivial zeros of the Riemann zeta function.

The program:
- loads Riemann gamma values from an embedded table or an external file;
- generates the corresponding prime sequence;
- solves the integer-lattice constraint
  a^2 + b^2 + c^2 + d^2 = p,
  a*d - b*c = D
  for each prime p and phase-derived starting value D;
- compares the Riemann-based reconstruction to a seeded random-phase baseline;
- exports a side-by-side CSV summary and Monte Carlo statistics;
- optionally runs a GUE semicircle unfolding check.

## What is being tested

The core idea is to map a prime-zero pair into a 4-square quaternion-style lattice representation, then assess whether the deterministic Riemann phase produces a tighter or more structured solution than a random phase baseline. The program measures exact matches, search steps, and Monte Carlo significance.

## Build

From this directory:

```bash
g++ -std=c++17 -O2 main.cpp -o spectral_inversion
```

## Run

```bash
./spectral_inversion --count 100 --max-steps 40 --mc-trials 500 --seed 42 --output spectral_comparison_100.csv
```

Common options:

```bash
./spectral_inversion --help
```

Key flags:
- `--gammas-file <path>`: load gamma values from a file instead of the built-in table
- `--count <N>`: number of prime-zero pairs to evaluate
- `--max-steps <N>`: maximum search radius for the lattice match
- `--mc-trials <N>`: Monte Carlo baseline trials
- `--seed <N>`: RNG seed
- `--output <path>`: CSV output path
- `--verify-gue`: run the GUE semicircle CDF unfolding check
- `--quiet`: suppress per-row terminal output

## Output

Running the program generates:
- a side-by-side CSV comparison file (default: `spectral_comparison_100.csv`);
- terminal statistics including exact matches, mean steps, Z-score, and Monte Carlo p-value;
- an optional GUE verification report when `--verify-gue` is enabled.

The repository already includes a sample run in `consoleOutput.txt`, which reflects a 100-node evaluation with default parameters.

## Files

- `main.cpp`: executable entry point and CLI driver
- `spectral_inversion.h`: shared constants, data structures, and default gamma table
- `lattice_solver.hpp`: prime generation and exact lattice solving
- `phase_accumulator.hpp`: split-phase accumulator and unit-circle validation
- `io_and_montecarlo.hpp`: CSV export and Monte Carlo comparison logic
- `gue_unfolding.hpp`: GUE semicircle unfolding check
- `rational_projection.hpp`: rational projection helpers
- `big_rational.hpp`: arbitrary-precision rational arithmetic support
- `spectral_comparison_100.csv`: generated comparison data for the default 100-node run
- `consoleOutput.txt`: sample console output from a representative execution

## Notes

- This is a numerical and symbolic prototype, not a production-ready library.
- The implementation is intentionally focused on exact integer arithmetic and explicit lattice verification.
- The built-in default gamma list captures the first several non-trivial zeros of the Riemann zeta function, and larger runs can be supplied through `--gammas-file`.
