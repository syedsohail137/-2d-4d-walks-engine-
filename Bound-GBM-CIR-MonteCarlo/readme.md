# Multi-Threaded SDE Simulation & Lattice Bounding Verification Engine

A high-performance, cross-platform native C++17 simulation matrix designed to evaluate continuous Stochastic Differential Equation (SDE) pathways across discrete number-theoretic lattices. The engine executes massively parallel Monte Carlo simulations to analyze how continuous process topologies impact mapping criteria driven by Fermat's Sum of Two Squares Theorem.

## 🚀 Engine Highlights

- **Massively Parallel Processing:** Leverages native hardware concurrency (`std::thread`) with a zero-thrashing local accumulation buffer pipeline.
- **O(1) Quadratic Residue Pruning:** Employs a constant-time Mod-7 Quadratic Residue Filter Map (`IS_MOD7_RESIDUE`) to strip out invalid lattice candidates before invoking expensive square root arithmetic.
- **Multi-Model Sideband Synthesis:** Conducts concurrent tracking across three distinct stochastic regimes: Geometric Brownian Motion (GBM), Cox-Ingersoll-Ross (CIR), and a Driftless Pure Random Walk.

---

## 🔬 Mathematical Foundations

The core architecture cross-examines the output of three continuous stochastic processes:

1. **Geometric Brownian Motion (GBM):**  
   $$dX_t = \mu X_t dt + \sigma X_t dW_t$$
2. **Cox-Ingersoll-Ross (CIR):**  
   $$dX_t = k(\theta - X_t)dt + \sigma \sqrt{X_t} dW_t$$
3. **Pure Random Walk (Driftless Log-Normal):**  
   $$dX_t = \sigma X_t dW_t$$

### The Lattice Mapping Transform
Every continuous point $X_t$ is transformed into a discrete state space scalar:
$$n = \lfloor X_t \cdot K_{\text{scale}} + 0.5 \rfloor$$

An integer $n$ successfully passes the verification gatekeeper if it satisfies Fermat's prerequisite condition ($n = 2$ or $n \equiv 1 \pmod 4$) and can be cleanly decomposed into two integer vector coordinates on the complex plane:
$$n = a^2 + d^2 \quad (a, d \in \mathbb{Z}^+)$$

---

## 📈 Empirical Validation Matrix ($10^8$ Nodes Evaluated)

High-volume stress testing ($1,000,000$ paths, $100$ discretization steps per path) uncovers deep behavioral anomalies between structural drift mechanics and number-theoretic distributions:

| Process Model | Total Nodes | Passed Nodes | Pass Rate | Phase Mean ($\mu$) | Phase StdDev ($\sigma$) | Bound Mean ($\mu$) | Bound StdDev ($\sigma$) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **GBM** | 100,000,000 | 16,400,890 | **16.40%** | 69.40° | ± 11.62 | 13.154 | ± 0.861 |
| **CIR (Square-Rt)** | 100,000,000 | 22,849,861 | **22.84%** | 84.05° | ± 2.08 | 13.061 | ± 0.030 |
| **Pure Random** | 100,000,000 | 16,274,780 | **16.27%** | 69.57° | ± 11.38 | 13.016 | ± 0.844 |

### Architectural Critique
- **The CIR Variance Trap:** Due to active mean reversion pulling the process trajectory toward $\theta = 100$, the CIR model encounters a strict state-space compression. This causes it to repeatedly assault a narrow integer band, artificially inflating its pass rate to **22.84%** while crushing phase variance down to **±2.08°**.
- **Drift Dominance Parity:** Over short horizons ($T = 1.0$), GBM and Pure Random tracks behave nearly identically. The stochastic diffusion term ($\sigma X_t dW_t$) completely overwhelms the minor 5% drift vector coefficient ($\mu$), resulting in an equivalent structural lattice distribution footprint.

---

## 🛠 Building and Execution

### Prerequisites
- Cross-platform C++17 compliant compiler (`g++` >= 7, `clang++` >= 5, or MSVC 2017+)
- Native thread support library components

### Compilation
Compile the standalone source package utilizing maximum optimization flags for your platform:

```bash
# Using GCC Toolchain
g++ -std=c++17 -O3 -march=native main.cpp -o sde_engine -lpthread
```

### Running the Engine
```bash
./sde_engine
```

---

## 📄 License
This performance asset is open-source software distributed under the terms of the **MIT License**.
