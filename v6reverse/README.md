# Riemann-Prime Spectral Lattice Projection Engine

A compact C++17 program that pairs the first non-trivial zeros of the Riemann zeta function with prime numbers and searches for integer 2x2 lattice coordinates that satisfy

$$
\begin{cases}
a^2+b^2+c^2+d^2 = p \\
ad-bc = D
\end{cases}
$$

for each prime $p$ and projected determinant $D$.

The implementation is intentionally lightweight: it generates primes, computes a phase angle from the Riemann zero, projects an initial determinant estimate, and searches nearby integer determinants until a valid lattice is found or the search limit is reached.

---

## Core idea

For the $n$-th non-trivial zero $\gamma_n$ and the $n$-th prime $p_n$, the engine computes a phase angle

$$
\theta_n = (\gamma_n \cdot \ln(p_n)) \bmod 2\pi
$$

and then uses a projected determinant

$$
D_{\text{initial}} = \mathrm{round}\left(\frac{p_n}{2}\cos(\theta_n)\right)
$$

as the starting point for the search.

The algebraic reduction is based on

$$
\begin{cases}
(a+d)^2 + (b-c)^2 = p + 2D \\
(a-d)^2 + (b+c)^2 = p - 2D
\end{cases}
$$

which lets the solver reduce the 4-variable problem to two sum-of-two-squares searches. The implementation finds all integer pairs for each target, expands sign/orientation combinations, and reconstructs candidate lattice coordinates.

---

## Codebase structure

- `spectral_inversion.h`: portable constants and data definitions, including `QuaternionLattice`, `SpectralNode`, and the `RIEMANN_GAMMAS` table for the first 30 non-trivial zeros.
- `main.cpp`: prime generation, phase projection, lattice solving, determinant search, invariant checking, and CSV export.
- `spectral_nodes.csv`: generated output containing the computed per-node state and lattice variables.

---

## Classification states

Every node is assigned one of the following states:

| State | Meaning |
| :--- | :--- |
| `Exact Match` | A valid lattice was found immediately at `D_initial`, with zero adjustment. |
| `Shifted Match` | No exact match existed at the initial determinant, so the solver searched outward to the nearest valid determinant. |
| `No Representation` | No valid lattice was found within the configured maximum radius. |

The program records:

- `Index`
- `Gamma_Full`
- `Prime`
- `Phase`
- `D_Initial`
- `D_Converged`
- `Search_Steps`
- `Adjustment`
- `State`
- `a, b, c, d`

---

## Build and run

### Prerequisites

Any C++17-compatible compiler such as `g++`, `clang++`, or MSVC.

### Build

```bash
g++ -O3 -std=c++17 -Wall -Wextra -pedantic main.cpp -o spectral_inversion
```

### Run

```bash
./spectral_inversion
```

This writes `spectral_nodes.csv` in the current working directory.

---

## Algorithm summary

1. Generate the first `n` primes using a sieve.
2. Load the first 30 non-trivial Riemann zeta zeros from `RIEMANN_GAMMAS`.
3. For each zero-prime pair, compute the phase angle and initial determinant estimate.
4. Solve
   
   $$a^2+b^2+c^2+d^2 = p$$
   
   and
   
   $$ad-bc = D$$
   
   using two sum-of-two-squares reductions.
5. If the initial determinant fails, search nearby integer values `D_initial ± 1, ±2, ...` within a bounded radius.
6. Verify the found lattice satisfies both the norm equation and determinant invariant before exporting CSV data.

---

## Current results

For the first 30 zero-prime pairs in the repository, the program resolves the dataset with the following behavior:

- `Exact Match`: 12 / 30
- `Shifted Match`: 18 / 30
- `No Representation`: 0 / 30

These values are produced by the current fixed table (`RIEMANN_GAMMAS`) and the project’s default `max_steps = 40` search window.

---

## Notes and scope

This is a proof-of-concept / research-oriented solver rather than a general-purpose theorem prover.

- It is built around a fixed set of the first 30 Riemann zeros.
- The determinant search is bounded and intentionally simple.
- The implementation is designed to validate the algebraic structure of the mapping, not to claim a universal spectral theorem.

The code is intentionally transparent: it does not hide the search logic, invariant checks, or CSV export path, and it exposes the exact determinant and lattice values it finds.

---

## References and interpretation

The implementation is based on a direct algebraic reduction from a 4-variable lattice problem to two quadratic forms.

The relevant mathematical structure is:

- prime numbers $p$
- non-trivial zeros $\gamma_n$ of the Riemann zeta function
- phase projection $\theta_n = \gamma_n \ln p_n \bmod 2\pi$
- determinant search around $D = \mathrm{round}(\frac{p}{2}\cos\theta)$
- lattice reconstruction from sum-of-two-squares representations

This project is best understood as an experimental numeric engine for testing structural correspondences between spectral quantities and integer lattice constraints.
