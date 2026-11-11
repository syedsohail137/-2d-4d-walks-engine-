# Riemann-Prime Spectral Lattice Projection Engine

This directory implements a compact C++17 prototype that pairs non-trivial Riemann zeta zeros with prime numbers and searches for integer 2x2 lattice coordinates satisfying

$$
\begin{cases}
 a^2 + b^2 + c^2 + d^2 = p \\
 ad - bc = D
\end{cases}
$$

where $p$ is the prime and $D$ is the projected determinant.

## What the code does

1. Generates the first $n$ primes using a simple sieve.
2. Loads the first 30 non-trivial Riemann zeta zeros from `RIEMANN_GAMMAS`.
3. Computes a phase angle from $\gamma_n \log p_n$.
4. Projects an initial determinant estimate.
5. Solves the lattice equations with a reduced two-squares formulation.
6. Searches nearby determinant values if an exact match is not found.
7. Writes CSV output containing the final state and lattice coefficients.

## Files

- `spectral_inversion.h` — constants and data structures
- `main.cpp` — lattice solver and CSV exporter
- `spectral_nodes.csv` — generated output from a run

## Build and run

```bash
g++ -O3 -std=c++17 -Wall -Wextra -pedantic main.cpp -o spectral_inversion
./spectral_inversion
```

This writes `spectral_nodes.csv` in the current working directory.

## State labels

The program records one of the following states:

- `Exact Match`
- `Shifted Match`
- `No Representation`

The output includes `Index`, `Gamma_Full`, `Prime`, `Phase`, `D_Initial`, `D_Converged`, `Search_Steps`, `Adjustment`, `State`, and the lattice values `(a, b, c, d)`.
