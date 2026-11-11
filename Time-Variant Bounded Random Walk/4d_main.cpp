#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

struct Config {
    int R = 8;                  // Number of local candidates to try per prime step
    double barrier = 1.0;       // Log-barrier strength; 0 disables barrier
    double bound = 50.0;        // Bounded walk radius; this is an empirical cap, not a theorem
    bool parity_twist = true;   // Apply (-1)^det sign flip when det is odd
    double phase_weight = 1.0;  // Weight of phase modulation; 0 disables it
    int N = 2000;               // Sieve limit
    int seed = 7;               // Kept for compatibility; not used by the current logic
};

struct Step {
    long long a = 0;
    long long b = 0;
    long long c = 0;
    long long d = 0;
    long long det = 0;
};

static bool four_square(long long p, Step& s) {
    if (p < 0) return false;
    long long limit = static_cast<long long>(std::sqrt(static_cast<double>(p)));
    for (long long a = -limit; a <= limit; ++a) {
        for (long long b = -limit; b <= limit; ++b) {
            for (long long c = -limit; c <= limit; ++c) {
                long long rem = p - a * a - b * b - c * c;
                if (rem < 0) continue;
                long long d = static_cast<long long>(std::sqrt(static_cast<double>(rem)));
                if (d * d != rem) continue;
                s = {a, b, c, d, a * d - b * c};
                return true;
            }
        }
    }
    return false;
}

static std::vector<Step> candidates(const Step& base, int R) {
    std::vector<Step> out;
    const long long v[4] = {base.a, base.b, base.c, base.d};

    static const int perms[24][4] = {
        {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
        {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
        {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
        {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}
    };

    for (int i = 0; i < 24 && static_cast<int>(out.size()) < R; ++i) {
        for (int mask = 0; mask < 16 && static_cast<int>(out.size()) < R; ++mask) {
            Step s{};
            long long vals[4] = {0, 0, 0, 0};
            for (int k = 0; k < 4; ++k) {
                vals[k] = v[perms[i][k]];
                if ((mask >> k) & 1) vals[k] = -vals[k];
            }
            s = {vals[0], vals[1], vals[2], vals[3], vals[0] * vals[3] - vals[1] * vals[2]};
            out.push_back(s);
        }
    }
    return out;
}

static double barrier_cost(double x, double y, const Config& cfg) {
    if (cfg.barrier == 0.0) return 0.0;

    double r = std::hypot(x, y);
    if (r >= cfg.bound) return 1e9;

    // As the radius approaches the empirical bound, the penalty grows rapidly.
    return -cfg.barrier * std::log(cfg.bound - r);
}

int main() {
    Config cfg;

    std::vector<int> is_composite(cfg.N + 1, 0);
    std::vector<long long> primes;
    primes.reserve(cfg.N / std::log(cfg.N));

    for (int i = 2; i <= cfg.N; ++i) {
        if (!is_composite[i]) {
            primes.push_back(i);
            for (long long j = 1LL * i * i; j <= cfg.N; j += i) {
                is_composite[static_cast<int>(j)] = 1;
            }
        }
    }

    const double gamma1 = 14.134725141734693;
    const double two_pi = 2.0 * std::acos(-1.0);

    double sx = 0.0;
    double sy = 0.0;

    std::printf("Step,Prime,Theta,Det,Parity,Sx,Sy,Norm,Barrier_Cost\n");

    int step = 0;
    for (long long p : primes) {
        Step base;
        if (!four_square(p, base)) continue;

        double theta = std::fmod(gamma1 * std::log(static_cast<double>(p)), two_pi);
        if (theta < 0.0) theta += two_pi;

        double px = base.a * std::cos(theta) * cfg.phase_weight +
                    base.b * std::sin(theta) * cfg.phase_weight +
                    base.a * (1.0 - cfg.phase_weight);
        double py = base.c * std::cos(theta) * cfg.phase_weight +
                    base.d * std::sin(theta) * cfg.phase_weight +
                    base.c * (1.0 - cfg.phase_weight);
        (void)px;
        (void)py;

        auto cands = candidates(base, cfg.R);
        double best_score = std::numeric_limits<double>::infinity();
        Step best = base;
        double best_x = sx;
        double best_y = sy;

        for (auto& c : cands) {
            long long det = c.det;
            long long sgn = 1;
            if (cfg.parity_twist && (det & 1LL)) sgn = -1;

            double nx = sx + sgn * static_cast<double>(c.a);
            double ny = sy + sgn * static_cast<double>(c.c);
            double norm = std::hypot(nx, ny);
            double score = norm + barrier_cost(nx, ny, cfg);

            if (score < best_score) {
                best_score = score;
                best = c;
                best_x = nx;
                best_y = ny;
            }
        }

        sx = best_x;
        sy = best_y;

        double norm = std::hypot(sx, sy);
        int parity = (cfg.parity_twist && (best.det & 1LL)) ? -1 : 1;

        std::printf("%d,%lld,%.6f,%lld,%d,%.3f,%.3f,%.3f,%.3f\n",
                    ++step,
                    p,
                    theta,
                    best.det,
                    parity,
                    sx,
                    sy,
                    norm,
                    barrier_cost(sx, sy, cfg));

        if (norm >= cfg.bound) {
            // Measurement cap reached. Do not claim proof beyond this radius.
        }
    }

    return 0;
}
