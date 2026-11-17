# Owen Gearbox: Symmetrized Scrambled Sobol + Adaptive Control Variate Prototype

This repository contains a compact C++17 numerical prototype for investigating a small set of variance-reduction ideas in a simple setting. In particular, it explores:

- a low-discrepancy Sobol-like sequence,
- a four-point symmetry transformation (`Sym4` / Gearbox), and
- an adaptive control-variate adjustment.

The project is intended as a research and teaching prototype rather than a validated production QMC library. It is useful for experimenting with ideas and for comparing a few implementation choices on a simple GBM/CIR benchmark, but it should not be interpreted as a general claim that the method is superior to standard Sobol/QMC approaches in all settings.

The current implementation includes several caveats that are worth keeping in mind:

- the Sobol-style generator here is a lower-triangular binary scramble, which is not automatically equivalent to a full Owen-type randomization unless explicitly verified;
- the CIR payoff is evaluated with a finite time discretization, so the estimated bond price includes both sampling and time-stepping error;
- the reported comparisons are based on a small number of fixed runs and should be read as illustrative examples, not as universal performance guarantees.

---

## 1. System architecture

The implementation separates the random-point generation step from the variance-reduction transformation step.

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

The code is intentionally compact and self-contained, with no external dependencies.

---

## 2. Mathematical setup

### 2.1 GBM payoff

For the standard GBM model,

- `S_T = S_0 exp((r - 0.5 sigma^2) T + sigma sqrt(T) Z)`,
- payoff = `exp(-r T) max(S_T - K, 0)`,
- control variate `C = S_T`,
- exact expectation = `E[C] = S_0 exp(r T)`.

This is the standard European call setup used here as a benchmark.

### 2.2 CIR bond payoff

For the CIR short-rate model,

- `dr_t = kappa (theta - r_t) dt + sigma sqrt(r_t) dW_t`,
- the exact transition can be represented using the standard Poisson-Gamma / noncentral-chi-square form,
- the code samples each time step using that transition law,
- the bond payoff is approximated by discounting over a finite time partition.

The benchmark uses the analytic bond price, while the simulation estimate still includes discretization error due to the finite time grid.

### 2.3 Adaptive control variate

For each batch,

- `Y_i = average of the transformed payoff values`,
- `C_i = average of the transformed control values`,
- `beta_hat = Cov(Y, C) / Var(C)`,
- adjusted estimator = `Ybar - beta_hat * (Cbar - E[C])`.

This is a standard control-variate estimator, but the coefficient is estimated from the same sample as the final estimate. That is a practical approximation and can introduce finite-sample bias, especially when the sample size is modest.

---

## 3. Sobol / Gearbox implementation notes

The code includes a scrambled Sobol-like sequence generator and a four-point symmetry transformation used to create a group orbit. This is a reasonable numerical experiment, but several caveats should be kept in mind:

1. The implementation is a prototype rather than a formally validated Sobol library.
2. The scramble is a lower-triangular binary linear scramble; it is not automatically the same as a standard Owen scrambling unless separately validated.
3. Results can depend strongly on the chosen randomization scheme and the number of independent replicates.
4. The Gearbox transformation is not guaranteed to improve every integrand; it may reduce variance for some structured problems but does not automatically dominate plain Sobol for all tasks.

This project should therefore be read as a compact exploratory method, not as a general result proving that the Gearbox construction is uniformly better than standard QMC methods.

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

The program prints the following benchmark summary for the fixed configuration used in the example:

```text
=== Owen Gearbox - GBM + CIR - Equal-Work Validation ===
GBM True BS   : 10.450584
CIR True Bond : 0.958791

GBM MC               | paths=4096 | est=10.440186 | sd=0.068094 | err=0.010397
GBM Pure Sobol       | paths=4096 | est=10.450848 | sd=0.002418 | err=0.000264
GBM Gearbox Sym4+CV  | paths=4096 | est=10.441864 | sd=0.000621 | err=0.008719

CIR MC               | paths=2048 | est=0.958835 | sd=0.000153 | err=0.000044
CIR Pure Sobol       | paths=2048 | est=0.958756 | sd=0.000024 | err=0.000035
CIR Gearbox Sym4+CV  | paths=2048 | est=0.958804 | sd=0.000005 | err=0.000013
```

This output is illustrative and should not be treated as a universal claim of superiority. In particular:

- for the GBM case, the Gearbox estimator has a smaller reported standard deviation but a larger absolute error than plain Sobol in this run;
- the code does not include a formal statistical validation criterion beyond a simple comparison of point estimates;
- the CIR result is encouraging, but it still depends on the finite time partition used in the simulation.

---

## 6. Recommended validation procedure

For a more rigorous and defensible QMC benchmark, the following checks are recommended:

1. Use a trusted randomized Sobol implementation or verify the direction numbers and randomization scheme against a reference source.
2. Include a proper digital shift or standard randomization step, and report independent replicates.
3. For CIR, test multiple time-step sizes and report convergence as the time grid is refined.
4. Use many randomized QMC replicates (for example, 32, 64, or more) and report confidence intervals for the estimator.
5. Estimate the control-variate coefficient independently or use out-of-sample estimation to reduce finite-sample bias.
6. Benchmark across multiple parameter sets and payoff functions rather than a single GBM and a single bond setup.

---

## 7. Summary

This project is a compact educational prototype for exploring Sobol-based variance reduction with a symmetry transformation and a control variate. The formulas and numerical ideas are reasonable, and the code is useful as a sandbox for experimentation and hypothesis testing.

At the same time, stronger validation is needed before claiming robust improvement over standard Sobol/QMC methods. The present implementation is best understood as a small exploratory benchmark rather than a definitive numerical result.
