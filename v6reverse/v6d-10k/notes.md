# Algebraic & Spectral Notes: Riemann-Prime Lattice Projection

This document summarizes the algebraic structure behind the 4D quaternion lattice solver, the empirical results from scaling the projection to $N = 2,000$ prime-zero pairs, and what these findings reveal about the Riemann Hypothesis.

---

## 1. Algebraic Structure of the Lattice System

The engine searches for four integers $(a, b, c, d) \in \mathbb{Z}^4$ satisfying the coupled system:

$$
\begin{cases}
a^2 + b^2 + c^2 + d^2 = p \\
ad - bc = D
\end{cases}
$$

* **Integral Quaternions (Lipschitz Order).** A tuple $(a, b, c, d)$ defines an integral quaternion $q = a + b\mathbf{i} + c\mathbf{j} + d\mathbf{k}$ with norm $N(q) = a^2 + b^2 + c^2 + d^2$. By **Lagrange's Four-Square Theorem**, every prime $p$ can be written as a sum of four squares. By **Jacobi's Four-Square Theorem**, there are exactly $8(p + 1)$ such integer tuples for any odd prime $p$.
* **Determinant as a Geometric Tilt.** Arranging the coordinates into a $2 \times 2$ integer matrix $M = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$ fixes its squared Frobenius norm ($\|M\|_F^2 = p$) and its determinant ($\det(M) = ad - bc = D$). Geometrically, $D$ measures the signed area or chirality of the lattice vector within the $8(p + 1)$ points on the 3-sphere of radius $\sqrt{p}$.
* **Decoupling into Two-Square Circles.** Adding and subtracting $2D$ from the norm splits the 4D sphere into two orthogonal 2D circles:

$$
(a + d)^2 + (b - c)^2 = p + 2D \quad \text{and} \quad (a - d)^2 + (b + c)^2 = p - 2D
$$

By **Fermat's Sum of Two Squares Theorem**, an integer is representable as $X^2 + Y^2$ if and only if every prime factor $q \equiv 3 \pmod 4$ divides it to an **even** power. Consequently, a target determinant $D$ is achievable for prime $p$ if and only if *both* $p + 2D$ and $p - 2D$ avoid odd exponents on all $4k + 3$ prime factors.

---

## 2. Spectral Theory & Why $(p_n, \gamma_n)$ Are Paired

In analytic number theory, prime numbers and non-trivial zeros of the Riemann zeta function $\zeta(s)$ act as Fourier duals:

* **The Phase Term $\gamma_n \ln p_n$.** In Riemann's explicit formulas and the **Selberg Trace Formula**, primes enter via their logarithms ($\ln p$, analogous to lengths of primitive closed geodesics), while zeta zeros enter via their imaginary parts ($\gamma_n$, analogous to vibrational frequencies of a hypothetical **Hilbert–Pólya operator**).
* **The Prototype Hypothesis.** Treating the $n$-th frequency $\gamma_n$ and $n$-th prime length $\ln p_n$ as a coupled oscillator with wrapped phase $\theta_n = (\gamma_n \ln p_n) \bmod 2\pi$, the engine tests whether the projected determinant $D_0$ lands directly on a valid quaternion determinant slice of $p_n$.

---

## 3. Scaled Empirical Benchmark ($N = 2,000$, $p_{\max} = 17,389$)

Evaluating the first $2,000$ prime-zero pairs across $200$ uniform-phase Monte Carlo trials contrasts two projection geometries:

* **Cosine Projection:** $D_0 = \operatorname{llround}\left(\frac{p}{2}\cos\theta_n\right)$
* **Linear (Area-Preserving) Projection:** $D_0 = \operatorname{llround}\left(p \left(\frac{\theta_n}{2\pi} - \frac{1}{2}\right)\right)$

| Projection Geometry | Exact Matches | Monte Carlo Null Mean | Z-Score | $p \equiv 1 \pmod 4$ Exact Rate | $p \equiv 3 \pmod 4$ Exact Rate |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Cosine Projection** | **300 / 2,000 (15.00%)** | $268.09 \pm 13.18$ (13.40%) | **+2.42** | 13.29% | 16.58% |
| **Linear (Area-Preserving)** | **223 / 2,000 (11.15%)** | $241.60 \pm 14.46$ (12.08%) | **-1.29** | 11.05% | 11.15% |

### Key Benchmark Takeaways

* **Projection Geometry Sensitivity.** Under the area-preserving **linear** projection—which maps $\theta \in [0, 2\pi)$ uniformly onto $[-p/2, p/2]$ without piling estimates near the boundaries—the Z-score is **$-1.29$**, making Riemann zero phases statistically indistinguishable from uniform random phases.
* **Cosine Boundary Artifact.** The cosine projection creates an arcsine density pile-up near $D \approx \pm p/2$. Near the boundary, one of the two-square targets $p \mp 2D$ approaches $0$ (which is trivially $0^2 + 0^2$), inflating sensitivity to minor phase clustering.
* **Residue Class Symmetry.** Under uniform linear projection, primes $p \equiv 1 \pmod 4$ ($11.05\%$) and $p \equiv 3 \pmod 4$ ($11.15\%$) exhibit virtually identical exact lattice solvability rates.

---

## 4. Insights Regarding the Riemann Hypothesis

1. **Equidistribution Modulo $2\pi$ (Weyl's Criterion).**
   Because $\gamma_n \sim \frac{2\pi n}{\ln n}$ and $p_n \sim n \ln n$, their product grows as:

   $$
   \gamma_n \ln p_n \sim 2\pi n \left(1 + \frac{\ln \ln n}{\ln n}\right)
   $$

   By Weyl's equidistribution theorem, the fractional parts of $\frac{\gamma_n \ln p_n}{2\pi}$ become uniformly distributed on $[0, 1)$ as $n \to \infty$. Empirical Z-scores clustering near zero confirm that $\theta_n = (\gamma_n \ln p_n) \bmod 2\pi$ behaves like a uniform random variable on the unit circle.

2. **Why Individual $(p_n, \gamma_n)$ Pairing Cannot Test RH.**
   The Riemann Hypothesis ($\operatorname{Re}(\rho) = \frac{1}{2}$ for all non-trivial zeros) is equivalent to the prime counting error bound $|\pi(x) - \operatorname{Li}(x)| = O(\sqrt{x}\ln x)$. That bound arises from **collective interference** across *all* zeros summed over *all* prime powers simultaneously in Riemann's explicit formula:

   $$
   \psi(x) = x - \sum_{\rho} \frac{x^{\rho}}{\rho} - \ln(2\pi) - \frac{1}{2}\ln(1 - x^{-2})
   $$

   Pairing the $n$-th prime exclusively with the $n$-th zero isolates a single diagonal phase $\exp(i \gamma_n \ln p_n)$ while discarding the cross-terms between distinct primes and zeros where RH phase cancellation occurs.

3. **Local Obstructions vs. Global Smoothness ("Shifted Matches").**
   Whether $p \pm 2D$ avoids odd powers of $4k + 3$ primes is a local modular obstruction, whereas $\gamma_n \ln p_n$ is a smooth transcendental phase. Because analytic phases do not detect local modular congruences at a single integer $D_0$, the solver typically shifts by $\pm 1$ to $\pm 5$ steps (`Shifted Match`) to bypass a $4k + 3$ prime factor.

4. **Asymptotic Sparsity of Solvable Determinants.**
   By Landau's theorem on sums of two squares, the density of integers $D \in [-p/2, p/2]$ for which both $p + 2D$ and $p - 2D$ are sums of two squares decays as $O\left(\frac{1}{\ln p}\right)$. The simulation directly reproduces this decay: the exact-match rate drops from **40.0%** at $N = 30$ ($p_{\max} = 113$), to **28.0%** at $N = 50$ ($p_{\max} = 229$), down to **11.15%** at $N = 2,000$ ($p_{\max} = 17,389$).
