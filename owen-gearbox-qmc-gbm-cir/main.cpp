// ============================================================
//  Owen Gearbox - GBM + CIR - Portable C++17 (Block 1 of 4)
//  Compile: g++ -O2 -std=c++17 main.cpp -o run && ./run
// ============================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace std;

// ---------- Acklam Inverse Normal CDF (|err| < 1.15e-9) ----------
static double norm_ppf(double p) {
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                             -2.759285104469687e+02, 1.383577518672690e+02,
                             -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                             -1.556989798598866e+02, 6.680131188771972e+01,
                             -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                             -2.400758277161838e+00, -2.549732539343734e+00,
                             4.374664141464968e+00,  2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                             2.445134137142996e+00, 3.754408661907416e+00};
  const double plow = 0.02425;
  const double phigh = 1.0 - plow;
  double q, r;
  if (p < plow) {
    q = sqrt(-2.0 * log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
           ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }
  if (p <= phigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
  }
  q = sqrt(-2.0 * log(1.0 - p));
  return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
         ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
}

// ---------- Incomplete Gamma & Noncentral Chi-Square PPF ----------
static double gammainc_reg(double a, double x) {
  if (x <= 0.0) return 0.0;
  if (x < a + 1.0) {
    double sum = 1.0 / a, del = sum;
    for (int n = 1; n < 1000; n++) {
      del *= x / (a + n);
      sum += del;
      if (fabs(del) < fabs(sum) * 1e-14) break;
    }
    return sum * exp(-x + a * log(x) - lgamma(a));
  }
  double b = x + 1.0 - a, c = 1e300, d = 1.0 / b, h = d;
  for (int i = 1; i < 1000; i++) {
    double an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs(d) < 1e-300) d = 1e-300;
    c = b + an / c;
    if (fabs(c) < 1e-300) c = 1e-300;
    d = 1.0 / d;
    double del = d * c;
    h *= del;
    if (fabs(del - 1.0) < 1e-14) break;
  }
  return 1.0 - exp(-x + a * log(x) - lgamma(a)) * h;
}

static double gamma_ppf(double u, double shape, double scale) {
  double lo = 0.0;
  double hi = shape * scale + 20.0 * sqrt(shape) * scale + 20.0;
  for (int i = 0; i < 64; i++) {
    double mid = 0.5 * (lo + hi);
    if (gammainc_reg(shape, mid / scale) < u)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

static double poisson_pmf(int k, double lam) {
  if (lam <= 0.0) return (k == 0) ? 1.0 : 0.0;
  return exp(-lam + k * log(lam) - lgamma(k + 1.0));
}

static int poisson_ppf(double u, double lam) {
  if (lam <= 0.0) return 0;
  double cdf = 0.0;
  for (int k = 0; k < 10000; k++) {
    cdf += poisson_pmf(k, lam);
    if (u <= cdf) return k;
  }
  return 10000;
}

static double ncx2_ppf(double u1, double u2, double dof, double nc) {
  int J = poisson_ppf(u1, nc * 0.5);
  double shape = (dof + 2.0 * J) * 0.5;
  return gamma_ppf(u2, shape, 2.0);
}
// ============================================================
//  Linear-Scrambled Sobol Sequence (Block 2 of 4)
// ============================================================
struct SobolDir {
  int s;
  int a;
  int m[18];
};

static const SobolDir DIRS[8] = {
    {1, 0, {1}},
    {2, 1, {1, 3}},
    {3, 1, {1, 3, 1}},
    {3, 2, {1, 3, 3}},
    {4, 1, {1, 3, 5, 13}},
    {4, 3, {1, 3, 5, 1}},
    {5, 2, {1, 3, 1, 3, 7}},
    {5, 4, {1, 3, 3, 9, 7}},
};

struct ScrambledSobol {
  int dim;
  vector<array<uint32_t, 32>> V;
  vector<array<array<int, 32>, 32>> L;

  ScrambledSobol(int dim_, uint64_t seed) : dim(dim_) {
    mt19937_64 sr(seed);
    V.assign(dim, {});
    for (int d = 0; d < dim; d++) {
      const auto& dir = DIRS[d];
      uint32_t m[32] = {0};
      for (int i = 0; i < dir.s; i++) m[i] = (uint32_t)dir.m[i] << (31 - i);
      for (int i = dir.s; i < 32; i++) {
        uint32_t v = m[i - dir.s] ^ (m[i - dir.s] >> dir.s);
        for (int k = 1; k < dir.s; k++) {
          v ^= (((dir.a >> (dir.s - 1 - k)) & 1) ? m[i - k] : 0);
        }
        m[i] = v;
      }
      for (int j = 0; j < 32; j++) V[d][j] = m[j];
    }
    L.assign(dim, {});
    for (int d = 0; d < dim; d++) {
      for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) {
          if (i == j)
            L[d][i][j] = 1;
          else if (j < i)
            L[d][i][j] = (sr() & 1);
          else
            L[d][i][j] = 0;
        }
      }
    }
  }

  double at(uint64_t n, int d) const {
    uint32_t x = 0;
    uint64_t g = n ^ (n >> 1);
    for (int j = 0; j < 32; j++) {
      if ((g >> j) & 1ULL) x ^= V[d][j];
    }
    uint32_t y = 0;
    for (int i = 0; i < 32; i++) {
      int bit = 0;
      for (int j = 0; j <= i; j++) bit ^= (L[d][i][j] & ((x >> (31 - j)) & 1));
      if (bit) y |= (1u << (31 - i));
    }
    // Midpoint offset guarantees symmetry under (1 - u) and prevents 0 or 1
    return ((y >> 8) + 0.5) * (1.0 / 16777216.0);
  }
};
// ============================================================
//  Analytical Truth, Adaptive CV, and GBM Engine (Block 3 of 4)
// ============================================================
static double bs_call(double S, double K, double r, double sig, double T) {
  double d1 = (log(S / K) + (r + 0.5 * sig * sig) * T) / (sig * sqrt(T));
  double d2 = d1 - sig * sqrt(T);
  auto N = [](double x) { return 0.5 * erfc(-x / sqrt(2.0)); };
  return S * N(d1) - K * exp(-r * T) * N(d2);
}

static double cir_bond(double r0, double kappa, double theta, double sig, double T) {
  double g = sqrt(kappa * kappa + 2.0 * sig * sig);
  double eg = exp(g * T);
  double B = 2.0 * (eg - 1.0) / ((g + kappa) * (eg - 1.0) + 2.0 * g);
  double A = pow(2.0 * g * exp((kappa + g) * T * 0.5) / ((g + kappa) * (eg - 1.0) + 2.0 * g),
                 2.0 * kappa * theta / (sig * sig));
  return A * exp(-B * r0);
}

static double apply_control_variate(const vector<double>& Y, const vector<double>& C, double EC) {
  size_t n = Y.size();
  double meanY = 0.0, meanC = 0.0;
  for (size_t i = 0; i < n; i++) {
    meanY += Y[i];
    meanC += C[i];
  }
  meanY /= n;
  meanC /= n;

  double cov = 0.0, varC = 0.0;
  for (size_t i = 0; i < n; i++) {
    double dc = C[i] - meanC;
    cov += (Y[i] - meanY) * dc;
    varC += dc * dc;
  }
  double beta = (varC > 1e-18) ? (cov / varC) : 0.0;
  return meanY - beta * (meanC - EC);
}

struct SimRes {
  double est;
  double sd;
  int evals;
};

static SimRes gbm(int total_paths, int R, int method) {
  const double S0 = 100.0, K = 100.0, r = 0.05, sig = 0.2, T = 1.0;
  const double EC = S0 * exp(r * T);
  vector<double> ests;
  int actual_evals = (method == 2) ? (total_paths / 4) * 4 : total_paths;

  for (int rep = 0; rep < R; rep++) {
    mt19937_64 mc_rng(42 + rep * 999);
    uniform_real_distribution<double> unif(1e-12, 1.0 - 1e-12);
    ScrambledSobol sob(1, 1234 + rep * 999);

    if (method != 2) {
      double sum = 0.0;
      for (int i = 0; i < actual_evals; i++) {
        double u = (method == 0) ? unif(mc_rng) : sob.at(i, 0);
        double Z = norm_ppf(u);
        double ST = S0 * exp((r - 0.5 * sig * sig) * T + sig * sqrt(T) * Z);
        sum += exp(-r * T) * max(ST - K, 0.0);
      }
      ests.push_back(sum / actual_evals);
    } else {
      int N_outer = actual_evals / 4;
      vector<double> Y_batch(N_outer), C_batch(N_outer);
      for (int i = 0; i < N_outer; i++) {
        double u = sob.at(i, 0);
        double u_shift = fmod(u + 0.5, 1.0);
        double us[4] = {u, 1.0 - u, u_shift, 1.0 - u_shift};
        double per = 0.0, C = 0.0;
        for (int k = 0; k < 4; k++) {
          double Z = norm_ppf(us[k]);
          double ST = S0 * exp((r - 0.5 * sig * sig) * T + sig * sqrt(T) * Z);
          per += exp(-r * T) * max(ST - K, 0.0);
          C += ST;
        }
        Y_batch[i] = per * 0.25;
        C_batch[i] = C * 0.25;
      }
      ests.push_back(apply_control_variate(Y_batch, C_batch, EC));
    }
  }

  double m = 0.0;
  for (double e : ests) m += e;
  m /= ests.size();
  double v = 0.0;
  for (double e : ests) v += (e - m) * (e - m);
  v /= max(1, (int)ests.size() - 1);
  return {m, sqrt(v), actual_evals};
}
// ============================================================
//  CIR Engine and Main Validation Harness (Block 4 of 4)
// ============================================================
static pair<double, double> cir_path(const vector<double>& U, double r0, double kappa,
                                     double theta, double sig, double T, int M) {
  double dt = T / M;
  double decay = exp(-kappa * dt);
  double c = sig * sig * (1.0 - decay) / (4.0 * kappa);
  double dof = 4.0 * kappa * theta / (sig * sig);

  double r = r0, disc = 0.0;
  for (int j = 0; j < M; j++) {
    double nc = 4.0 * kappa * decay * r / (sig * sig * (1.0 - decay));
    double rn = c * ncx2_ppf(U[2 * j], U[2 * j + 1], dof, nc);
    disc += 0.5 * (r + rn) * dt;
    r = rn;
  }
  return {exp(-disc), r};
}

static SimRes cir(int total_paths, int R, int method) {
  const double kappa = 0.5, theta = 0.05, sig = 0.1, r0 = 0.04, T = 1.0;
  const int M = 4, DIM = 2 * M;
  const double EC = theta + (r0 - theta) * exp(-kappa * T);
  vector<double> ests;
  int actual_evals = (method == 2) ? (total_paths / 4) * 4 : total_paths;

  for (int rep = 0; rep < R; rep++) {
    mt19937_64 mc_rng(777 + rep * 131);
    uniform_real_distribution<double> unif(1e-12, 1.0 - 1e-12);
    ScrambledSobol sob(DIM, 777 + rep * 131);

    if (method != 2) {
      double sum = 0.0;
      vector<double> U(DIM);
      for (int i = 0; i < actual_evals; i++) {
        for (int d = 0; d < DIM; d++) U[d] = (method == 0) ? unif(mc_rng) : sob.at(i, d);
        sum += cir_path(U, r0, kappa, theta, sig, T, M).first;
      }
      ests.push_back(sum / actual_evals);
    } else {
      int N_outer = actual_evals / 4;
      vector<double> Y_batch(N_outer), C_batch(N_outer);
      vector<double> U(DIM), Uk(DIM);

      for (int i = 0; i < N_outer; i++) {
        for (int d = 0; d < DIM; d++) U[d] = sob.at(i, d);
        double per = 0.0, C = 0.0;
        for (int k = 0; k < 4; k++) {
          for (int d = 0; d < DIM; d++) {
            double u = U[d];
            double u_shift = fmod(u + 0.5, 1.0);
            if (k == 0)
              Uk[d] = u;
            else if (k == 1)
              Uk[d] = 1.0 - u;
            else if (k == 2)
              Uk[d] = u_shift;
            else
              Uk[d] = 1.0 - u_shift;
          }
          auto [payoff, rT] = cir_path(Uk, r0, kappa, theta, sig, T, M);
          per += payoff;
          C += rT;
        }
        Y_batch[i] = per * 0.25;
        C_batch[i] = C * 0.25;
      }
      ests.push_back(apply_control_variate(Y_batch, C_batch, EC));
    }
  }

  double m = 0.0;
  for (double e : ests) m += e;
  m /= ests.size();
  double v = 0.0;
  for (double e : ests) v += (e - m) * (e - m);
  v /= max(1, (int)ests.size() - 1);
  return {m, sqrt(v), actual_evals};
}

int main() {
  cout.setf(ios::fixed);
  cout << setprecision(6);

  double true_bs = bs_call(100.0, 100.0, 0.05, 0.2, 1.0);
  double true_cir = cir_bond(0.04, 0.5, 0.05, 0.1, 1.0);

  cout << "=== Owen Gearbox - GBM + CIR - Equal-Work Validation ===\n";
  cout << "GBM True BS   : " << true_bs << "\n";
  cout << "CIR True Bond : " << true_cir << "\n\n";

  for (auto [name, method] :
       vector<pair<string, int>>{{"MC", 0}, {"Pure Sobol", 1}, {"Gearbox Sym4+CV", 2}}) {
    auto g = gbm(4096, 8, method);
    double err = fabs(g.est - true_bs);
    cout << "GBM " << left << setw(16) << name << " | paths=" << g.evals << " | est=" << g.est
         << " | sd=" << g.sd << " | err=" << err
         << (method == 2 && err < 0.01 ? "  PASS" : "") << "\n";
  }

  cout << "\n";

  for (auto [name, method] :
       vector<pair<string, int>>{{"MC", 0}, {"Pure Sobol", 1}, {"Gearbox Sym4+CV", 2}}) {
    auto c = cir(2048, 8, method);
    double err = fabs(c.est - true_cir);
    cout << "CIR " << left << setw(16) << name << " | paths=" << c.evals << " | est=" << c.est
         << " | sd=" << c.sd << " | err=" << err
         << (method == 2 && err < 0.005 ? "  PASS" : "") << "\n";
  }

  return 0;
}
