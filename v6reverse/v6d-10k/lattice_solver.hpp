#ifndef LATTICE_SOLVER_HPP
#define LATTICE_SOLVER_HPP

#include "spectral_inversion.h"

#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>

namespace spectral {

// Dynamic Sieve of Eratosthenes guaranteed to produce at least n primes
inline std::vector<long long> generate_primes(std::size_t n) {
    if (n == 0) return {};
    long long limit = (n < 6)
        ? 15
        : static_cast<long long>(n * (std::log(static_cast<double>(n)) +
                                      std::log(std::log(static_cast<double>(n))))) + 64;

    while (true) {
        std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
        std::vector<long long> primes;
        primes.reserve(n);

        for (long long p = 2; p <= limit; ++p) {
            if (is_prime[static_cast<std::size_t>(p)]) {
                primes.push_back(p);
                if (primes.size() == n) return primes;
                if (p <= limit / p) {
                    for (long long i = p * p; i <= limit; i += p) {
                        is_prime[static_cast<std::size_t>(i)] = false;
                    }
                }
            }
        }
        limit *= 2;
    }
}

// Finds all non-negative integer pairs (X, Y) such that X^2 + Y^2 = target
inline void solve_two_squares_into(long long target,
                                   std::vector<std::pair<long long, long long>>& out_pairs) {
    out_pairs.clear();
    if (target < 0) return;

    long long left = 0;
    long long right = static_cast<long long>(std::sqrt(static_cast<double>(target)));

    while (left <= right) {
        long long sum = left * left + right * right;
        if (sum == target) {
            out_pairs.emplace_back(right, left);
            if (right != left) {
                out_pairs.emplace_back(left, right);
            }
            ++left;
            --right;
        } else if (sum < target) {
            ++left;
        } else {
            --right;
        }
    }
}

// Thread-safe, reusable-buffer O(sqrt(p)) Diophantine solver
inline bool solve_lattice_buffered(long long p, long long D, QuaternionLattice& out,
                                   std::vector<std::pair<long long, long long>>& pairs1,
                                   std::vector<std::pair<long long, long long>>& pairs2) {
    if (std::llabs(D) > p / 2) return false;

    solve_two_squares_into(p + 2 * D, pairs1);
    if (pairs1.empty()) return false;

    solve_two_squares_into(p - 2 * D, pairs2);
    if (pairs2.empty()) return false;

    for (const auto& p1 : pairs1) {
        long long X1 = p1.first;
        long long Y1 = p1.second;
        for (const auto& p2 : pairs2) {
            long long X2 = p2.first;
            long long Y2 = p2.second;

            for (long long sx1 : {1LL, -1LL}) {
                if (X1 == 0 && sx1 == -1) continue;
                for (long long sy1 : {1LL, -1LL}) {
                    if (Y1 == 0 && sy1 == -1) continue;
                    for (long long sx2 : {1LL, -1LL}) {
                        if (X2 == 0 && sx2 == -1) continue;
                        for (long long sy2 : {1LL, -1LL}) {
                            if (Y2 == 0 && sy2 == -1) continue;

                            long long x1 = sx1 * X1;
                            long long y1 = sy1 * Y1;
                            long long x2 = sx2 * X2;
                            long long y2 = sy2 * Y2;

                            if (((x1 + x2) % 2 == 0) && ((y1 + y2) % 2 == 0)) {
                                long long a = (x1 + x2) / 2;
                                long long d = (x1 - x2) / 2;
                                long long b = (y1 + y2) / 2;
                                long long c = (y2 - y1) / 2;

                                if (a * a + b * b + c * c + d * d == p &&
                                    a * d - b * c == D) {
                                    out = {a, b, c, d};
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

inline SpectralNode resolve_from_initial_d(std::size_t idx, double gamma, long long p,
                                           double theta, long long d_initial,
                                           long long max_steps,
                                           std::vector<std::pair<long long, long long>>& buf1,
                                           std::vector<std::pair<long long, long long>>& buf2) {
    SpectralNode node;
    node.index = static_cast<int>(idx + 1);
    node.riemann_zero_gamma = gamma;
    node.prime = p;
    node.phase_angle_rad = theta;
    node.d_initial = d_initial;

    if (solve_lattice_buffered(p, d_initial, node.q, buf1, buf2)) {
        node.d_converged = d_initial;
        node.search_steps = 0;
        node.adjustment = 0;
        node.state = "Exact Match";
        return node;
    }

    for (long long step = 1; step <= max_steps; ++step) {
        for (long long delta : {step, -step}) {
            long long candidate_d = d_initial + delta;
            if (std::llabs(candidate_d) > p / 2) continue;

            if (solve_lattice_buffered(p, candidate_d, node.q, buf1, buf2)) {
                node.d_converged = candidate_d;
                node.search_steps = static_cast<int>(step);
                node.adjustment = std::llabs(candidate_d - d_initial);
                node.state = "Shifted Match";
                return node;
            }
        }
    }

    node.d_converged = d_initial;
    node.search_steps = static_cast<int>(max_steps);
    node.adjustment = -1;
    node.state = "No Representation";
    node.q = {0, 0, 0, 0};
    return node;
}

} // namespace spectral

#endif // LATTICE_SOLVER_HPP
