# Riemann-Prime Spectral Lattice Projection Engine

## Included Files

- `spectral_inversion.h` — Data structures, projection geometry helpers (`cosine` and `linear`), and embedded first 50 non-trivial Riemann zeta zeros.
- `lattice_solver.hpp` — Dynamic Sieve of Eratosthenes and zero-allocation $O(\sqrt{p})$ Diophantine quaternion lattice solver.
- `io_and_montecarlo.hpp` — External LMFDB/Odlyzko zero file parser, OpenMP-parallel Monte Carlo significance engine, and side-by-side CSV exporter.
- `main.cpp` — CLI entrypoint, algebraic invariant verifier, magnitude scale-tier breakdown, and timing diagnostics.
- `plot_spectral_results.py` — 4-panel Matplotlib diagnostic visualization script.

## Build & Run

```bash
# Compile with OpenMP support
g++ -O3 -std=c++17 -fopenmp -march=native -Wall -Wextra -pedantic main.cpp -o spectral_inversion

# Run default 50 embedded zeros with 1,000 Monte Carlo trials
./spectral_inversion

# Fetch Odlyzko's first 100,000 Riemann zeros and run 10,000-node sweep
curl -sL https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1 -o zeros1.txt
./spectral_inversion --gammas-file zeros1.txt --count 10000 --projection cosine --mc-trials 1000 --quiet --output run_cosine_10k.csv
./spectral_inversion --gammas-file zeros1.txt --count 10000 --projection linear --mc-trials 1000 --quiet --output run_linear_10k.csv

# Generate diagnostic plots
python3 plot_spectral_results.py run_cosine_10k.csv
python3 plot_spectral_results.py run_linear_10k.csv
```
