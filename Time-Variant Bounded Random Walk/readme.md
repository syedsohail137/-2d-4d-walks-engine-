# Time-Variant Bounded Random Walk

This directory contains experiments for bounded walk dynamics over prime-generated or number-theoretic coordinate spaces. The emphasis is on reducing discrepancy and stabilizing motion under structural constraints.

## Main ideas

- 2D and 4D bounded walk variants
- prime-lattice coordinate generation
- discrepancy minimization and barrier-like penalty logic
- deterministic or structured update rules instead of purely unconstrained random motion

## Build

Prerequisites:
- C++17 or later
- OpenMP-enabled compiler recommended
- CMake optional if used in a broader build flow

Compile:

```bash
g++ -O3 -std=c++17 -fopenmp main.cpp -o bounded_walk_engine
```

Run:

```bash
./bounded_walk_engine
```

## Notes

This directory contains experimental code and related data outputs; it is meant for exploratory research rather than as a production optimizer.
