# One Gearbox: Quasi-Monte Carlo Variance Reduction Engine

A zero-dependency C++17 benchmark demonstrating combined variance reduction for numerical integration and option pricing. It pairs Owen-scrambled Sobol sequences with a 4-point symmetry orbit (`Sym4`) and exact analytical control variates.

## Methodology

The engine evaluates each sample path across three stages:

1. **Owen-Scrambled Sobol Sequence:** Generates low-discrepancy points in $[0, 1)^D$ using 24-dimensional Joe-Kuo direction numbers (`new-joe-kuo-6.21201`) with hash-based bit scrambling across $R = 16$ independent replicates ($N = 4096$ base points).
2. **`Sym4` Orbit:** Maps each uniform coordinate $u$ to the 4-element orbit $\{u, 1-u, (u+0.5) \bmod 1, 1 - ((u+0.5) \bmod 1)\}$ before applying Acklam's inverse normal CDF approximation.
3. **Exact Control Variate:** Adjusts the orbit-averaged payoff $f$ using a correlated control $C$ with a known closed-form expectation $\mathbb{E}[C]$ via ordinary least squares:
   $$\hat{\mu} = \frac{1}{N} \sum_{i=1}^{N} \left( \bar{f}_i - \hat{\beta}(\bar{C}_i - \mathbb{E}[C]) \right)$$

## Benchmark Results ($N = 4096, R = 16$)

All models run with $S_0 = 100$, $K = 100$, $r = 0.05$, $\sigma = 0.2$, $T = 1.0$, and $12$ time steps.

| Benchmark | Control Variate ($C$) | Truth | Plain MC SD | Gearbox SD | VRF | \|Bias\| |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **GBM Call** | Terminal price $S_T$ | `10.4506` | `0.2060` | `0.0006` | **138,565x** | `0.0001` |
| **CIR Bond** | Terminal short rate $r_T$ | `0.9513` | `0.0002` | `0.0000` | **225x** | `0.0000` |
| **Asian Arithmetic** | Geometric Asian payoff | `6.1516` | `0.1288` | `0.0019` | **4,597x** | `0.0043` |
| **Heston Call** | Terminal price $S_T$ | `10.4298` | `0.1563` | `0.0756` | **4.3x** | `0.0305` |
| **American Put (LSM)** | None (`Sym4` + Owen only) | `6.0310` | `0.1403` | `0.0506` | **7.7x** | `0.0638` |
| **Kepler Mean $E$** | None (`Sym4` exact) | `3.1416` | — | — | — | `0.0000` |

*Note: VRF (Variance Reduction Factor) is computed as $(\text{SD}_{\text{MC}} / \text{SD}_{\text{OURS}})^2$ across 16 replicates.*

## Known Limitations

* **Time-Discretization Bias (Heston & CIR):** Full-truncation Euler stepping with 12 monthly steps leaves a residual discretization bias relative to continuous-time analytical truths (`0.0305` on Heston).
* **In-Sample LSM Bias (American Put):** Least-Squares Monte Carlo regression is fitted and evaluated on the same sample paths without a control variate, introducing a small upward exercise bias (`0.0638`).
* **Dimensionality:** Direction numbers are included for the first 24 dimensions ($D \le 24$). Higher-dimensional problems require extending `JK_TABLE`.

## Build & Run

Requires a standard C++17 compiler with no external libraries.

```bash
g++ -O3 -std=c++17 -Wall -Wextra main.cpp -o gearbox
./gearbox
```
## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

MIT License

Copyright (c) 2026 [Full Name]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

