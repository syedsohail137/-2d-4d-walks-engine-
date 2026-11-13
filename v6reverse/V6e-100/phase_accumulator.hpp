#ifndef PHASE_ACCUMULATOR_HPP
#define PHASE_ACCUMULATOR_HPP

#include "lattice_solver.hpp"
#include "spectral_inversion.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <utility>

namespace spectral {

struct PrimeSplit {
    long long Ra = 0, Rb = 0, Rc = 0, Rd = 0;
    long long n1 = 0, n2 = 0;
    long long det = 0;
};

struct PhaseAccumulatorState {
    double zr = 1.0;
    double zi = 0.0;
    double total = 0.0;
    BigRational q_re{BigInt(1), BigInt(1)};
    BigRational q_im{BigInt(0), BigInt(1)};
};

inline PrimeSplit best_prime_split(long long p) {
    PrimeSplit best{};
    long long best_abs_det = -1;

    for (long long n1 = 1; n1 < p; ++n1) {
        long long n2 = p - n1;
        auto v1 = solve_two_squares_all(n1);
        auto v2 = solve_two_squares_all(n2);
        if (v1.empty() || v2.empty()) continue;

        auto expand = [](const std::vector<std::pair<long long, long long>>& raw) {
            std::vector<std::pair<long long, long long>> out;
            for (const auto& pr : raw) {
                out.push_back(pr);
                if (pr.first != pr.second) out.emplace_back(pr.second, pr.first);
            }
            return out;
        };

        auto p1 = expand(v1);
        auto p2 = expand(v2);

        for (auto [a, b] : p1) {
            for (auto [c, d] : p2) {
                for (long long sb : {1LL, -1LL}) {
                    for (long long sd : {1LL, -1LL}) {
                        long long det = a * (sd * d) - (sb * b) * c;
                        if (std::llabs(det) > best_abs_det) {
                            best_abs_det = std::llabs(det);
                            best = {a, sb * b, c, sd * d, n1, n2, det};
                        }
                    }
                }
            }
        }
    }
    return best;
}

inline long long sign_fix(long long Rb, long long p) {
    return (p % 4 == 3) ? -Rb : Rb;
}

inline std::tuple<double, double, double> normalize_pair(double a, double b) {
    double n = a * a + b * b;
    if (n == 0.0) return {1.0, 0.0, 0.0};
    double s = std::sqrt(n);
    return {a / s, b / s, std::atan2(b, a)};
}

// Exact S^1(Q) rational squaring lift: ((a^2 - b^2)/(a^2 + b^2), 2ab/(a^2 + b^2))
inline std::pair<BigRational, BigRational> rational_unit_lift_big(long long a, long long b) {
    long long n1 = a * a + b * b;
    if (n1 == 0) return {BigRational(BigInt(1), BigInt(1)), BigRational(BigInt(0), BigInt(1))};
    return {
        BigRational(BigInt(a * a - b * b), BigInt(n1)),
        BigRational(BigInt(2 * a * b), BigInt(n1))
    };
}

inline double wrap_to_pi(double angle) {
    double w = std::fmod(angle + PI, 2.0 * PI);
    if (w < 0.0) w += 2.0 * PI;
    return w - PI;
}

inline SplitPhaseRecord advance_split_phase(
    PhaseAccumulatorState& st,
    long long p,
    double riemann_phase_rad) {

    PrimeSplit sp = best_prime_split(p);
    long long b_fixed = sign_fix(sp.Rb, p);

    auto [zr, zi, theta] = normalize_pair(static_cast<double>(sp.Ra), static_cast<double>(b_fixed));

    // Brahmagupta-Fibonacci multiplication on S^1
    double new_re = st.zr * zr - st.zi * zi;
    double new_im = st.zr * zi + st.zi * zr;
    double norm = std::hypot(new_re, new_im);
    new_re /= norm;
    new_im /= norm;

    // Arbitrary-precision Brahmagupta-Fibonacci multiplication on S^1(Q)
    auto [ur, ui] = rational_unit_lift_big(sp.Ra, b_fixed);
    BigRational next_q_re = (st.q_re * ur) - (st.q_im * ui);
    BigRational next_q_im = (st.q_re * ui) + (st.q_im * ur);

    double total2 = st.total + theta;
    double wrapped = wrap_to_pi(total2);

    st.zr = new_re;
    st.zi = new_im;
    st.total = total2;
    st.q_re = next_q_re;
    st.q_im = next_q_im;

    SplitPhaseRecord rec;
    rec.Ra = sp.Ra;
    rec.Rb = sp.Rb;
    rec.Rc = sp.Rc;
    rec.Rd = sp.Rd;
    rec.split_det = sp.det;
    rec.step_theta = theta;
    rec.cumulative_total = total2;
    rec.cumulative_wrapped = wrapped;
    rec.zr = new_re;
    rec.zi = new_im;
    rec.q_re = next_q_re;
    rec.q_im = next_q_im;
    rec.riemann_phase_diff = wrap_to_pi(riemann_phase_rad - wrapped);
    return rec;
}

inline bool verify_split_phase_invariants(long long p, const SplitPhaseRecord& rec) {
    long long norm4 = rec.Ra * rec.Ra + rec.Rb * rec.Rb + rec.Rc * rec.Rc + rec.Rd * rec.Rd;
    if (norm4 != p) {
        std::cerr << "Split 4-square invariant failure at p=" << p << "\n";
        return false;
    }
    if (std::abs(std::hypot(rec.zr, rec.zi) - 1.0) > 1e-12) {
        std::cerr << "S^1 unit circle drift failure at p=" << p << "\n";
        return false;
    }
    BigRational norm_sq = (rec.q_re * rec.q_re) + (rec.q_im * rec.q_im);
    if (!norm_sq.num.is_one() || !norm_sq.den.is_one()) {
        std::cerr << "BigRational S^1(Q) unit norm failure at p=" << p << "\n";
        return false;
    }
    return true;
}

} // namespace spectral

#endif // PHASE_ACCUMULATOR_HPP
