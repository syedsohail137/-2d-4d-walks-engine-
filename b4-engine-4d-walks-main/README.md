# Project Blueprint: 4D Random Walk Engine (`B4_Engine`)

This directory contains one of the primary 4D walk prototypes in the repository. It is a C++ development project for experimenting with random walks, scaling behavior, and prime-based state generation in a higher-dimensional lattice setting.

## Overview

The project explores:

- four-dimensional random-walk dynamics
- coefficient and boundary handling
- prime-based or arithmetic-constrained state generation
- empirical drift and stability observations

## Structure

- `src/` — engine implementation files
- `include/` — headers and shared declarations
- `v2/` — earlier variant
- `V3/` — later variant

## Build

From this directory:

```bash
g++ -O3 -std=c++17 src/main.cpp -Iinclude -o b4_engine
```

## Run

```bash
./b4_engine
```

## Notes

This is an experimental codebase. It is useful for exploring ideas and comparing variants, but it should not be assumed to be a general-purpose library or a fully stabilized production engine.
