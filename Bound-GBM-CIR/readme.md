# Multi-Threaded SDE Simulation & Lattice Bounding Verification Engine

A high-performance, parallelized C++17 application designed to simulate Stochastic Differential Equations (SDEs) and filter continuous paths through a discrete arithmetic lattice.

## Features

* **Dual SDE Models:** Parallel simulation of Geometric Brownian Motion (GBM) and Cox-Ingersoll-Ross (CIR) processes.
* **Lattice Mapping:** Discretizes continuous paths and verifies them using Fermat’s Sum of Two Squares theorem (n = a² + d²).
* **O(1) Mod-7 Pre-Filtering:** Employs a quadratic residue filter to discard invalid coordinates instantly before expensive square root computations.
* **Native Threading:** Scalable multi-threaded execution utilizing local thread buffers to minimize global synchronization overhead.
* **Statistical Reporting:** Outputs unified metrics detailing total/passed nodes, pass rates, complex-plane phase distributions, and spectral tracking bounds.

---

## Technical Specifications

### SDE Configurations
* **GBM:** Multiplicative noise drift with a safety trap at zero.
* **CIR:** Mean-reverting drift featuring a square-root asset volatility structure and reflecting boundary enforcement.

### Verification Criteria
* **Fermat Condition:** n = 2 or \(n \equiv 1 \pmod 4\).
* **Vector Phase (°):** \(\theta = \arctan2(d, a) \times \frac{180}{\pi}\).
* **Spectral Tracking Bound:** Bound \(= 1.3 \times \sqrt{n}\).

---

## Build and Run

### Prerequisites
* A standard C++17 or higher compatible compiler (`g++`, `clang++`, or `MSVC`).

### Compilation
```bash
g++ -O3 -std=c++17 main.cpp -o sde_simulation -pthread
```

### Execution
```bash
./sde_simulation
```

---

## Example Performance Profile

```text
================================================================================
               SIDE-BY-SIDE SIDEBAND STRUCTURAL REPORT MATRIX                   
================================================================================
Process Model   Total Nodes    Passed Nodes   Pass Rate    Phase (μ ± σ)   Bound (μ ± σ)
------------------------------------------------------------------------------------------
GBM             1000000        163973         16.39%       69.44° ± 11.60  13.152 ± 0.859
CIR (Square-Rt) 1000000        227794         22.77%       84.05° ± 2.08   13.061 ± 0.030
------------------------------------------------------------------------------------------
Execution Status: 100% SDE Integration Verification Success.
================================================================================
```

---

## Architecture Design Notes

* **`SDEParameters`**: Encapsulates variables for initial valuation ($X_0$), drift coefficients ($\mu$, $k$), standard deviations ($\sigma$), bounds ($T$), and discrete scale coefficients ($K_{\text{scale}}$).
* **`AggregateStats`**: Structure storing accumulated moments (sum and sum of squares) to guarantee accurate standard deviations across millions of simulation nodes without keeping complete path histories in memory.
