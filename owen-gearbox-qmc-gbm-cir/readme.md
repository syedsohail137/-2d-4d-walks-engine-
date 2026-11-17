# Owen Gearbox: Symmetrized Scrambled Sobol + Adaptive Control Variate Engine

This repository contains a compact C++17 numerical prototype for exploring the interaction between:

- a low-discrepancy Sobol sequence,
- a four-point symmetry transformation (`Sym4` / Gearbox), and
- an adaptive control-variate adjustment.

It is intended as a research/demo implementation, not as a production-validated QMC library.

The current implementation is useful for studying variance reduction ideas in a simple setting, but some claims in earlier drafts should be read carefully. In particular:

- the Sobol randomization here is a lower-triangular linear scramble, not necessarily a full modern Owen-type randomization unless separately verified;
- the CIR payoff is evaluated with a finite time discretization, so the reported bond estimate contains both sampling and time-stepping error;
- the reported validation is based on a small number of fixed replicates and is best read as an illustrative benchmark rather than proof of universal superiority.

---

## 1. System architecture

The engine separates the generation of random-like points from the variance-reduction transformations applied to them.

```text
+---------------------------------------------------------------+
| 1. LOW-DISCREPANCY POINT GENERATOR                             |
| Sobol-type sequence with a binary scrambling step               |
|                                                               |
|  -> Generates points u in (0, 1)^d                            |
+---------------------------------------------------------------+
                              |
                              v
+---------------------------------------------------------------+
| 2. GEARBOX TRANSFORMATION                                      |
| Sym4-style group orbit:                                         |
|   u^(1) = u                                                   |
|   u^(2) = 1 - u                                               |
|   u^(3) = (u + 0.5) mod 1                                     |
|   u^(4) = 1 - ((u + 0.5) mod 1)                               |
+---------------------------------------------------------------+
                              |
                              v
+---------------------------------------------------------------+
| 3. PATH KERNELS                                               |
| GBM European call payoff                                       |
| CIR bond payoff                                                |
+---------------------------------------------------------------+
                              |
                              v
+---------------------------------------------------------------+
| 4. ADAPTIVE CONTROL VARIATE                                    |
| beta_hat = Cov(Y, C) / Var(C)                                  |
| final estimate = Ybar - beta_hat * (Cbar - E[C])              |
+---------------------------------------------------------------+
```

This implementation is intentionally compact and self-contained, with zero external dependencies.

---

## 2. Mathematical setup

### 2.1 GBM payoff

For the standard GBM model,

- S_T = S_0 exp((r - 0.5 sigma^2) T + sigma sqrt(T) Z),
- payoff = exp(-r T) max(S_T - K, 0),
- control variate C = S_T,
- exact expectation = E[C] = S_0 exp(r T).

This is the standard European call setup and is used for the benchmark.

### 2.2 CIR bond payoff

For the CIR short-rate model,

- dr_t = kappa (theta - r_t) dt + sigma sqrt(r_t) dW_t,
- the exact transition can be represented in the standard Poisson-Gamma / noncentral-chi-square form,
- the code samples each step using that transition law,
- the bond payoff is approximated by discounting over a finite time partition.

The exact analytic bond price is used as the benchmark, but the simulation estimate is still affected by the finite discretization in time.

### 2.3 Adaptive control variate

For each batch,

- Y_i = average of the transformed payoff values,
- C_i = average of the transformed control values,
- beta_hat = Cov(Y, C) / Var(C),
- adjusted estimator = Ybar - beta_hat * (Cbar - E[C]).

This is a standard control-variate estimator, but the coefficient is estimated from the same sample used to form the final estimate, which is a practical approximation rather than a formal finite-sample proof.

---

## 3. Sobol / Gearbox implementation notes

The code includes a scrambled Sobol-like sequence generator and a four-point transformation used to create a symmetry group orbit. This is a reasonable numerical experiment, but several caveats are important:

1. The implementation should be regarded as a prototype, not as a formally verified, production Sobol library.
2. The scramble is a lower-triangular binary linear scramble; this is not automatically equivalent to a standard Owen scrambling unless separately validated.
3. The reported variation across replicates depends strongly on the chosen randomization scheme and the number of independent replicates.
4. The Gearbox transformation is not guaranteed to improve every integrand; it can reduce variance in structured problems but does not automatically dominate plain Sobol on all tasks.

The README should therefore be read as describing a specific exploratory method, not as a general theorem about QMC superiority.

---

## 4. Build and run

```bash
g++ -O2 -std=c++17 main.cpp -o run && ./run
```

Requirements:

- any standard C++17 compiler,
- no external libraries required.

---

## 5. Example output from the current implementation

The program prints the following benchmark summary for the given fixed setup:

```text
=== Owen Gearbox - GBM + CIR - Equal-Work Validation ===
GBM True BS   : 10.450584
CIR True Bond : 0.958791

GBM MC               | paths=4096 | est=10.440186 | sd=0.068094 | err=0.010397
GBM Pure Sobol       | paths=4096 | est=10.450848 | sd=0.002418 | err=0.000264
GBM Gearbox Sym4+CV  | paths=4096 | est=10.441864 | sd=0.000621 | err=0.008719  PASS

CIR MC               | paths=2048 | est=0.958835 | sd=0.000153 | err=0.000044
CIR Pure Sobol       | paths=2048 | est=0.958756 | sd=0.000024 | err=0.000035
CIR Gearbox Sym4+CV  | paths=2048 | est=0.958804 | sd=0.000005 | err=0.000013  PASS
```

This output is illustrative and should not be treated as a universal claim of superiority. In particular:

- for the GBM case, the Gearbox estimator has a smaller reported standard deviation but a larger absolute error than pure Sobol in this run;
- the code marks a result as `PASS` using a simple heuristic threshold (`err < 0.01`), which is not a rigorous validation criterion;
- the CIR case is much more encouraging, but it still includes discretization error from the finite time partition.

---

## 6. Recommended validation procedure

For a more rigorous and defensible QMC benchmark, the following checks are recommended:

1. Use a trusted randomized Sobol implementation or verify the direction numbers and randomization scheme against a reference source.
2. Include a proper digital shift or standard randomization step, and report independent replicates.
3. For CIR, test multiple time-step sizes and report convergence as the time grid is refined.
4. Use many randomized QMC replicates (e.g. 32, 64, or more) and report confidence intervals for the estimator.
5. Estimate the control-variate coefficient independently or use out-of-sample estimation to reduce finite-sample bias.
6. Benchmark across multiple parameter sets and payoff functions, not just a single option and a single bond setup.

---

## 7. Summary

This project is a compact educational prototype for exploring Sobol-based variance reduction with a symmetry transformation and a control variate. The formulas and numerical ideas are reasonable, and the example outputs are plausible, but the current implementation and presentation should be interpreted with caution. The code is best viewed as a demonstration of a method under test rather than a formally validated pricing engine.

This is a useful sandbox for experimentation and hypothesis testing, but stronger validation is needed before claiming robust improvement over plain Sobol/QMC methods.
