# Discrete Lyapunov Optimization

This directory contains a control-oriented optimization experiment built around a discrete-state Lyapunov-style feedback mechanism.

## Goal

The project explores whether a structured transformation or optimization rule can reduce uncontrolled drift in a high-dimensional discrete system while preserving stability and boundedness.

## Typical implementation details

- multi-threaded C++ engine
- optimization over structured transformation matrices
- comparison between uncontrolled and controlled trajectories
- CSV/log output for convergence and stability diagnostics

## Build

Linux/macOS:

```bash
g++ -O3 -fopenmp -march=native main.cpp -o lyapunov_control_engine
```

Windows (MSVC):

```cmd
cl /O2 /openmp /arch:AVX2 main.cpp /Fe:lyapunov_control_engine.exe
```

## Run

```bash
./lyapunov_control_engine
```

## Notes

This directory is a standalone numerical experiment and is not intended as a general-purpose control library.
