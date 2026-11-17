#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

constexpr double PI = 3.14159265358979323846;

// ---------- Fast Xorshift64 PRNG ----------
static uint64_t rng_state = 0x123456789abcdefULL;
static inline uint64_t xorshift64() {
  uint64_t x = rng_state;
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  rng_state = x;
  return x;
}
static inline double rand01() {
  return (xorshift64() >> 11) * (1.0 / 9007199254740992.0);
}

// ---------- Fixed Normal PPF (Acklam) & CDF ----------
static double norm_ppf(double p) {
  static const double a[] = {
      -3.969683028665376e+01,  2.209460984245205e+02,
      -2.759285104469687e+02,  1.383577518672690e+02,
      -3.066479806614716e+01,  2.506628277459239e+00};
  static const double b[] = {
      -5.447609879822406e+01,  1.615858368580409e+02,
      -1.556989798598866e+02,  6.680131188771972e+01,
      -1.328068155288572e+01};
  static const double c[] = {
      -7.784894002430293e-03, -3.223964580411365e-01,
      -2.400758277161838e+00, -2.549732539343734e+00,
       4.374664141464968e+00,  2.938163982698783e+00};
  static const double d[] = {
       7.784695709041462e-03,  3.224671290700398e-01,
       2.445134137142996e+00,  3.754408661907416e+00};

  constexpr double plow = 0.02425;
  constexpr double phigh = 1.0 - plow;

  if (p < plow) {
    double q = sqrt(-2.0 * log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
           ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }
  if (p > phigh) {
    double q = sqrt(-2.0 * log(1.0 - p));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }
  double q = p - 0.5;
  double r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
         (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}

static inline double norm_cdf(double x) {
  return 0.5 * erfc(-x / sqrt(2.0));
}
// ---------- Full Joe-Kuo Sobol (24 dims, 32 bits) ----------
struct JK {
  int s;
  int a;
  uint32_t m[8];
};

static const JK JK_TABLE[24] = {
    {1, 0, {1}},
    {2, 1, {1, 3}},
    {3, 1, {1, 3, 7}},
    {3, 2, {1, 1, 5}},
    {4, 1, {1, 3, 1, 1}},
    {4, 6, {1, 1, 3, 7}},
    {5, 1, {1, 3, 3, 9, 9}},
    {5, 4, {1, 3, 7, 13, 3}},
    {5, 7, {1, 1, 5, 11, 27}},
    {5, 11, {1, 1, 7, 3, 29}},
    {5, 13, {1, 3, 7, 7, 21}},
    {5, 14, {1, 1, 1, 9, 23}},
    {5, 19, {1, 3, 1, 15, 17}},
    {5, 21, {1, 1, 5, 5, 1}},
    {5, 22, {1, 3, 3, 3, 25}},
    {5, 25, {1, 1, 3, 15, 29}},
    {5, 26, {1, 3, 1, 7, 3}},
    {5, 28, {1, 3, 7, 9, 31}},
    {5, 31, {1, 3, 1, 9, 5}},
    {6, 1, {1, 1, 1, 1, 1, 33}},
    {6, 3, {1, 1, 5, 3, 15, 19}},
    {6, 7, {1, 3, 1, 1, 9, 27}},
    {6, 11, {1, 1, 3, 7, 31, 15}},
    {6, 13, {1, 1, 3, 7, 3, 23}},
};

static uint32_t JK_V[24][32];
static bool JK_INIT = false;

static void jk_init() {
  if (JK_INIT) return;
  JK_INIT = true;
  for (int d = 0; d < 24; d++) {
    int s = JK_TABLE[d].s;
    int a = JK_TABLE[d].a;
    uint32_t m[32] = {0};
    for (int k = 0; k < s; k++) {
      m[k] = JK_TABLE[d].m[k];
    }
    for (int k = s; k < 32; k++) {
      // Fixed: include (m[k - s] << s) per Sobol recurrence
      uint32_t val = (m[k - s] << s) ^ m[k - s];
      for (int i = 1; i < s; i++) {
        int ai = (a >> (s - 1 - i)) & 1;
        if (ai) {
          val ^= (m[k - i] << i);
        }
      }
      m[k] = val;
    }
    for (int k = 0; k < 32; k++) {
      JK_V[d][k] = m[k] << (31 - k);
    }
  }
}

static inline uint32_t owen_scramble_bits(uint32_t x, uint64_t seed, int dim) {
  uint32_t out = 0;
  uint64_t base = seed * 0x9e3779b97f4a7c15ULL + (uint64_t)(dim + 1) * 1000003ULL;
  for (int j = 0; j < 32; j++) {
    int bit = (x >> (31 - j)) & 1;
    uint64_t h = base + (uint64_t)(out) * 131ULL + (uint64_t)j;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    if (h & 1) bit ^= 1;
    out = (out << 1) | (uint32_t)bit;
  }
  return out;
}

static inline uint32_t owen_sobol_uint32(int dim, int n, uint64_t seed) {
  uint32_t x = 0;
  int k = 0;
  int nn = n;
  while (nn) {
    if (nn & 1) x ^= JK_V[dim % 24][k];
    nn >>= 1;
    k++;
  }
  return owen_scramble_bits(x, seed + (uint64_t)dim * 7919ULL, dim);
}
// ---------- Simulation Constants & Closed-Form Benchmarks ----------
constexpr int N = 4096;
constexpr int R = 16;
constexpr int Msteps = 12;
constexpr double S0 = 100.0, K = 100.0, r = 0.05, sig = 0.2, T = 1.0;

static double bs_call(double S, double K_, double r_, double s_, double T_) {
  double d1 = (log(S / K_) + (r_ + 0.5 * s_ * s_) * T_) / (s_ * sqrt(T_));
  double d2 = d1 - s_ * sqrt(T_);
  return S * norm_cdf(d1) - K_ * exp(-r_ * T_) * norm_cdf(d2);
}

static double bs_put(double S, double K_, double r_, double s_, double T_) {
  return bs_call(S, K_, r_, s_, T_) - S + K_ * exp(-r_ * T_);
}

static double geo_asian_price() {
  double dt = T / Msteps;
  double sig2 = 0.0, mu = 0.0;
  for (int i = 1; i <= Msteps; i++) {
    double w = (double)(Msteps - i + 1) / Msteps;
    mu += w;
    sig2 += w * w;
  }
  double m = log(S0) + (r - 0.5 * sig * sig) * dt * mu;
  double v = sig * sig * dt * sig2;
  double sG = sqrt(v);
  double d1 = (m - log(K) + v) / sG;
  double d2 = d1 - sG;
  return exp(-r * T) * (exp(m + 0.5 * v) * norm_cdf(d1) - K * norm_cdf(d2));
}

static double cir_EC_rT() {
  double kappa = 0.5, theta = 0.05, r0 = 0.05;
  return theta + (r0 - theta) * exp(-kappa * T);
}

static double cir_bond_true() {
  double kappa = 0.5, theta = 0.05, s = 0.1, r0 = 0.05;
  double g = sqrt(kappa * kappa + 2.0 * s * s);
  double B = 2.0 * (exp(g * T) - 1.0) / ((g + kappa) * (exp(g * T) - 1.0) + 2.0 * g);
  double A = pow(2.0 * g * exp((kappa + g) * T / 2.0) /
                     ((g + kappa) * (exp(g * T) - 1.0) + 2.0 * g),
                 2.0 * kappa * theta / (s * s));
  return A * exp(-B * r0);
}

// ---------- Payoffs on Flattened Z Buffer (N * D) ----------
static void payoff_gbm(const vector<double>& Z, int D, vector<double>& f, vector<double>& C) {
  double disc = exp(-r * T);
  double drift = (r - 0.5 * sig * sig) * T;
  double vol = sig * sqrt(T);
  for (int i = 0; i < N; i++) {
    double ST = S0 * exp(drift + vol * Z[i * D]);
    f[i] = disc * max(ST - K, 0.0);
    C[i] = ST;
  }
}

static void payoff_asian(const vector<double>& Z, int D, vector<double>& f, vector<double>& C) {
  double dt = T / Msteps;
  double disc = exp(-r * T);
  double drift = (r - 0.5 * sig * sig) * dt;
  double vol = sig * sqrt(dt);
  for (int i = 0; i < N; i++) {
    double S = S0, A = 0.0, logsum = 0.0;
    const double* row = &Z[i * D];
    for (int j = 0; j < Msteps; j++) {
      S *= exp(drift + vol * row[j]);
      A += S;
      logsum += log(S);
    }
    A /= Msteps;
    f[i] = disc * max(A - K, 0.0);
    double G = exp(logsum / Msteps);
    C[i] = disc * max(G - K, 0.0);
  }
}

static void payoff_heston(const vector<double>& Z, int D, vector<double>& f, vector<double>& C) {
  double dt = T / Msteps, sqrt_dt = sqrt(dt);
  double disc = exp(-r * T);
  double v0 = 0.04, kappa = 2.0, theta = 0.04, xi = 0.3, rho = -0.7;
  double rho_comp = sqrt(1.0 - rho * rho);
  for (int i = 0; i < N; i++) {
    double S = S0, v = v0;
    const double* z1 = &Z[i * D];
    const double* z2 = &Z[i * D + Msteps];
    for (int j = 0; j < Msteps; j++) {
      double vv = max(v, 0.0);
      double sq_v = sqrt(vv);
      double dW1 = z1[j] * sqrt_dt;
      double dW2 = (rho * z1[j] + rho_comp * z2[j]) * sqrt_dt;
      S *= exp((r - 0.5 * vv) * dt + sq_v * dW1);
      v += kappa * (theta - v) * dt + xi * sq_v * dW2;
      v = max(v, 0.0);
    }
    f[i] = disc * max(S - K, 0.0);
    C[i] = S;
  }
}

static void payoff_cir(const vector<double>& Z, int D, vector<double>& f, vector<double>& C) {
  double kappa = 0.5, theta = 0.05, s = 0.1, r0 = 0.05;
  double dt = T / Msteps, sqrt_dt = sqrt(dt);
  for (int i = 0; i < N; i++) {
    double rr = r0, integ = 0.0;
    const double* row = &Z[i * D];
    for (int j = 0; j < Msteps; j++) {
      rr = max(rr, 0.0);
      integ += rr * dt;
      rr += kappa * (theta - rr) * dt + s * sqrt(rr) * sqrt_dt * row[j];
    }
    f[i] = exp(-integ);
    C[i] = max(rr, 0.0);
  }
}

static double solve_ls(const vector<double>& X, const vector<double>& Y, double x) {
  double s0 = (double)X.size(), s1 = 0, s2 = 0, s3 = 0, s4 = 0, t0 = 0, t1 = 0, t2 = 0;
  for (size_t i = 0; i < X.size(); i++) {
    double xi = X[i], yi = Y[i];
    double x2 = xi * xi;
    s1 += xi; s2 += x2; s3 += x2 * xi; s4 += x2 * x2;
    t0 += yi; t1 += xi * yi; t2 += x2 * yi;
  }
  double A[3][4] = {{s0, s1, s2, t0}, {s1, s2, s3, t1}, {s2, s3, s4, t2}};
  for (int col = 0; col < 3; col++) {
    int piv = col;
    for (int rr = col + 1; rr < 3; rr++) {
      if (fabs(A[rr][col]) > fabs(A[piv][col])) piv = rr;
    }
    swap(A[col], A[piv]);
    double d = A[col][col];
    if (fabs(d) < 1e-12) return 0.0;
    for (int cc = col; cc < 4; cc++) A[col][cc] /= d;
    for (int rr = 0; rr < 3; rr++) {
      if (rr != col) {
        double fac = A[rr][col];
        for (int cc = col; cc < 4; cc++) A[rr][cc] -= fac * A[col][cc];
      }
    }
  }
  return A[0][3] + A[1][3] * x + A[2][3] * x * x;
}

static void payoff_american(const vector<double>& Z, int D, vector<double>& f, vector<double>& C) {
  double dt = T / Msteps;
  double disc = exp(-r * dt);
  double drift = (r - 0.5 * sig * sig) * dt;
  double vol = sig * sqrt(dt);

  static vector<double> Sp(N * (Msteps + 1));
  for (int i = 0; i < N; i++) {
    Sp[i * (Msteps + 1)] = S0;
    const double* zrow = &Z[i * D];
    double* srow = &Sp[i * (Msteps + 1)];
    for (int j = 0; j < Msteps; j++) {
      srow[j + 1] = srow[j] * exp(drift + vol * zrow[j]);
    }
  }

  vector<double> V(N);
  double discT = exp(-r * T);
  for (int i = 0; i < N; i++) {
    double ST = Sp[i * (Msteps + 1) + Msteps];
    V[i] = max(K - ST, 0.0);
    C[i] = discT * V[i];
  }

  vector<double> X, Y;
  X.reserve(N);
  Y.reserve(N);

  for (int t = Msteps - 1; t >= 0; t--) {
    X.clear();
    Y.clear();
    for (int i = 0; i < N; i++) {
      V[i] *= disc;
      double st = Sp[i * (Msteps + 1) + t];
      if (K - st > 0.0) {
        X.push_back(st);
        Y.push_back(V[i]);
      }
    }
    if (X.size() > 20) {
      for (int i = 0; i < N; i++) {
        double st = Sp[i * (Msteps + 1) + t];
        double exercise = K - st;
        if (exercise > 0.0) {
          double cont = solve_ls(X, Y, st);
          if (exercise > cont) V[i] = exercise;
        }
      }
    }
  }
  f = V;
}
// ---------- Statistics & Gearbox Runners ----------
struct Res {
  double mean, sd, se, bias;
};

static Res stats(const vector<double>& v, double truth) {
  double m = 0.0;
  for (double x : v) m += x;
  m /= (double)v.size();
  double var = 0.0;
  for (double x : v) var += (x - m) * (x - m);
  var /= max(1.0, (double)v.size() - 1.0);
  double sd = sqrt(var), se = sd / sqrt((double)v.size());
  return {m, sd, se, fabs(m - truth)};
}

template <typename Payoff>
static vector<double> gearbox_run(Payoff payoff, int D, double EC, bool useBF) {
  jk_init();
  vector<double> rep_est(R);
  vector<double> U_base(N * D);
  vector<double> Z(N * D);
  vector<double> f(N), C(N), per(N), cmean(N);

  for (int rep = 0; rep < R; rep++) {
    uint64_t seed = 1000ULL + (uint64_t)rep * 7919ULL;
    fill(per.begin(), per.end(), 0.0);
    fill(cmean.begin(), cmean.end(), 0.0);

    // Precompute Owen-scrambled base points once per replicate
    for (int i = 0; i < N; i++) {
      for (int d = 0; d < D; d++) {
        uint32_t y = owen_sobol_uint32(d, i, seed);
        U_base[i * D + d] = (y + 0.5) / 4294967296.0;
      }
    }

    // Apply Sym4 orbit transforms on cached U_base
    for (int orb = 0; orb < 4; orb++) {
      for (int idx = 0; idx < N * D; idx++) {
        double u = U_base[idx];
        if (orb == 1) {
          u = 1.0 - u;
        } else if (orb == 2) {
          u = u + 0.5;
          if (u >= 1.0) u -= 1.0;
        } else if (orb == 3) {
          u = u + 0.5;
          if (u >= 1.0) u -= 1.0;
          u = 1.0 - u;
        }
        u = min(max(u, 1e-12), 1.0 - 1e-12);
        Z[idx] = norm_ppf(u);
      }

      payoff(Z, D, f, C);
      for (int i = 0; i < N; i++) {
        per[i] += f[i] * 0.25;
        cmean[i] += C[i] * 0.25;
      }
    }

    double beta = 0.0;
    if (useBF) {
      double mf = 0.0, mc = 0.0;
      for (int i = 0; i < N; i++) {
        mf += per[i];
        mc += cmean[i];
      }
      mf /= N;
      mc /= N;
      double cov = 0.0, var = 0.0;
      for (int i = 0; i < N; i++) {
        double df = per[i] - mf;
        double dc = cmean[i] - mc;
        cov += df * dc;
        var += dc * dc;
      }
      beta = (var > 1e-18 ? cov / var : 0.0);
    }

    double est = 0.0;
    for (int i = 0; i < N; i++) {
      est += per[i] - beta * (cmean[i] - EC);
    }
    rep_est[rep] = est / N;
  }
  return rep_est;
}

template <typename Payoff>
static vector<double> mc_run(Payoff payoff, int D) {
  vector<double> rep_est(R);
  vector<double> Z(N * D);
  vector<double> f(N), C(N);

  for (int rep = 0; rep < R; rep++) {
    rng_state = 0x123456789abcdefULL + (uint64_t)rep * 99991ULL;
    for (int idx = 0; idx < N * D; idx++) {
      double u = min(max(rand01(), 1e-12), 1.0 - 1e-12);
      Z[idx] = norm_ppf(u);
    }
    payoff(Z, D, f, C);
    double m = 0.0;
    for (double x : f) m += x;
    rep_est[rep] = m / N;
  }
  return rep_est;
}

int main() {
  cout.setf(ios::fixed);
  cout << setprecision(4);

  double truth_gbm = bs_call(S0, K, r, sig, T);
  double truth_cir = cir_bond_true();
  double truth_geo = geo_asian_price();
  double truth_asian = 6.1516, truth_heston = 10.4298, truth_amer = 6.0310;

  auto run_and_print = [&](const string& name, const vector<double>& mc,
                           const vector<double>& ours, double truth) {
    Res a = stats(mc, truth), b = stats(ours, truth);
    cout << name << " Truth=" << truth << "\n";
    cout << "  MC   mean=" << a.mean << " SD=" << a.sd << " SE=" << a.se
         << " |Bias|=" << a.bias << "\n";
    cout << "  OURS mean=" << b.mean << " SD=" << b.sd << " SE=" << b.se
         << " |Bias|=" << b.bias;
    if (a.sd > 0.0 && b.sd > 0.0) {
      cout << " VRF=" << (a.sd * a.sd) / (b.sd * b.sd) << "x";
    }
    cout << "\n";
    return b;
  };

  Res r_gbm = run_and_print(
      "GBM Call", mc_run(payoff_gbm, 1),
      gearbox_run(payoff_gbm, 1, S0 * exp(r * T), true), truth_gbm);

  Res r_cir = run_and_print(
      "CIR Bond", mc_run(payoff_cir, Msteps),
      gearbox_run(payoff_cir, Msteps, cir_EC_rT(), true), truth_cir);

  Res r_asian = run_and_print(
      "Asian Arith", mc_run(payoff_asian, Msteps),
      gearbox_run(payoff_asian, Msteps, truth_geo, true), truth_asian);

  Res r_heston = run_and_print(
      "Heston Call", mc_run(payoff_heston, 24),
      gearbox_run(payoff_heston, 24, S0 * exp(r * T), true), truth_heston);

  Res r_amer = run_and_print(
      "American Put (LSM)", mc_run(payoff_american, Msteps),
      gearbox_run(payoff_american, Msteps, bs_put(S0, K, r, sig, T), false),
      truth_amer);

  // Kepler Mean E over 4096 points with Sym4 orbit
  {
    int NK = 4096;
    double e = 0.5, sum = 0.0;
    for (int i = 0; i < NK; i++) {
      double U = (i + 0.5) / NK;
      double u_shift = (U + 0.5 >= 1.0) ? (U - 0.5) : (U + 0.5);
      double orbit[4] = {U, 1.0 - U, u_shift, 1.0 - u_shift};
      double s = 0.0;
      for (int k = 0; k < 4; k++) {
        double M = 2.0 * PI * orbit[k];
        double E = M + e * sin(M);
        E += (M - E + e * sin(E)) / (1.0 - e * cos(E));
        s += E;
      }
      sum += s * 0.25;
    }
    double meanE = sum / NK;
    cout << "Kepler Mean E Truth=" << PI << " Ours=" << meanE
         << " |Bias|=" << fabs(meanE - PI) << "\n";
  }

  bool pass = (r_gbm.bias < 0.01) && (r_cir.bias < 0.01) && (r_asian.bias < 0.02);
  cout << (pass ? "PASS\n" : "FAIL\n");
  return pass ? 0 : 1;
}
