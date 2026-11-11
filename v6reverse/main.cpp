#include "spectral_inversion.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

namespace {

// Sieve of Eratosthenes generating the first n primes
std::vector<long long> generate_primes(std::size_t n) {
    if (n == 0) return {};
    long long limit = (n < 6)
        ? 15
        : static_cast<long long>(n * (std::log(static_cast<double>(n)) +
                                      std::log(std::log(static_cast<double>(n)))));

    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    std::vector<long long> primes;
    primes.reserve(n);

    for (long long p = 2; p < limit + 1; ++p) {
        if (is_prime[static_cast<std::size_t>(p)]) {
            primes.push_back(p);
            if (primes.size() == n) break;
            for (long long i = p * p; i < limit + 1; i += p) {
                is_prime[static_cast<std::size_t>(i)] = false;
            }
        }
    }
    return primes;
}

// Finds all non-negative integer pairs (X, Y) with X >= Y such that X^2 + Y^2 = target
std::vector<std::pair<long long, long long>> solve_two_squares_all(long long target) {
    if (target < 0) return {};
    std::vector<std::pair<long long, long long>> solutions;
    long long left = 0;
    long long right = static_cast<long long>(std::sqrt(static_cast<double>(target)));

    while (left < right + 1) {
        long long sum = left * left + right * right;
        if (sum == target) {
            solutions.emplace_back(right, left);
            ++left;
            --right;
        } else if (sum < target) {
            ++left;
        } else {
            --right;
        }
    }
    return solutions;
}

// Exact O(sqrt(p)) Diophantine solver for:
//   a^2 + b^2 + c^2 + d^2 = p
//   a*d - b*c = D
bool solve_lattice(long long p, long long D, spectral::QuaternionLattice& out) {
    // Cauchy-Schwarz bound: |D| <= p / 2
    if (std::llabs(D) > p / 2) return false;

    auto raw1 = solve_two_squares_all(p + 2 * D);
    auto raw2 = solve_two_squares_all(p - 2 * D);
    if (raw1.empty() || raw2.empty()) return false;

    // Expand (X, Y) to include (Y, X) so all parity alignments are tested
    auto expand_pairs = [](const std::vector<std::pair<long long, long long>>& raw) {
        std::vector<std::pair<long long, long long>> pairs;
        pairs.reserve(raw.size() * 2);
        for (const auto& pr : raw) {
            pairs.push_back(pr);
            if (pr.first != pr.second) {
                pairs.emplace_back(pr.second, pr.first);
            }
        }
        return pairs;
    };

    const auto pairs1 = expand_pairs(raw1);
    const auto pairs2 = expand_pairs(raw2);

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

spectral::SpectralNode invert_node(std::size_t idx, long long p,
                                   long long max_steps, bool use_random_baseline) {
    spectral::SpectralNode node;
    node.index = static_cast<int>(idx + 1);
    node.riemann_zero_gamma = spectral::RIEMANN_GAMMAS[idx];
    node.prime = p;

    double theta = std::fmod(node.riemann_zero_gamma * std::log(static_cast<double>(p)),
                             2.0 * spectral::PI);
    if (theta < 0.0) theta += 2.0 * spectral::PI;
    node.phase_angle_rad = theta;

    if (use_random_baseline) {
        unsigned long long hash =
            static_cast<unsigned long long>(idx) * 2654435761ULL +
            static_cast<unsigned long long>(p);
        long long span = p / 2;
        node.d_initial = static_cast<long long>(hash % static_cast<unsigned long long>(2 * span + 1)) - span;
    } else {
        node.d_initial = std::llround(static_cast<double>(p) * 0.5 * std::cos(theta));
    }

    // Step 1: Direct determinant projection
    if (solve_lattice(p, node.d_initial, node.q)) {
        node.d_converged = node.d_initial;
        node.search_steps = 0;
        node.adjustment = 0;
        node.state = "Exact Match";
        return node;
    }

    // Step 2: Bounded nearest-integer determinant search
    for (long long step = 1; step < max_steps + 1; ++step) {
        for (long long delta : {step, -step}) {
            long long candidate_d = node.d_initial + delta;
            if (std::llabs(candidate_d) > p / 2) continue;

            if (solve_lattice(p, candidate_d, node.q)) {
                node.d_converged = candidate_d;
                node.search_steps = static_cast<int>(step);
                node.adjustment = std::llabs(candidate_d - node.d_initial);
                node.state = "Shifted Match";
                return node;
            }
        }
    }

    // Step 3: No representation found within max_steps
    node.d_converged = node.d_initial;
    node.search_steps = static_cast<int>(max_steps);
    node.adjustment = -1;
    node.state = "No Representation";
    node.q = {0, 0, 0, 0};
    return node;
}

} // namespace

int main() {
    const std::size_t count = spectral::RIEMANN_GAMMAS.size();
    const std::vector<long long> primes = generate_primes(count);
    const long long max_steps = 40;

    std::vector<spectral::SpectralNode> nodes;
    nodes.reserve(count);

    for (std::size_t i = 0; i < count; ++i) {
        nodes.push_back(invert_node(i, primes[i], max_steps, false));
    }

    std::ofstream out_csv("spectral_nodes.csv");
    if (!out_csv) {
        std::cerr << "Error: Failed to open spectral_nodes.csv for writing.\n";
        return 1;
    }

    const std::string header =
        "Index,Gamma_Full,Prime,Phase,D_Initial,D_Converged,Search_Steps,Adjustment,State,a,b,c,d\n";
    out_csv << header;
    std::cout << header;

    for (const auto& n : nodes) {
        // Verify algebraic invariants before output
        if (n.state != "No Representation") {
            long long norm_sq = n.q.a * n.q.a + n.q.b * n.q.b + n.q.c * n.q.c + n.q.d * n.q.d;
            long long det = n.q.a * n.q.d - n.q.b * n.q.c;
            if (norm_sq != n.prime || det != n.d_converged) {
                std::cerr << "Verification failure at Index " << n.index << "\n";
                return 2;
            }
        }

        out_csv << std::fixed
                << n.index << ","
                << std::setprecision(15) << n.riemann_zero_gamma << ","
                << n.prime << ","
                << std::setprecision(6) << n.phase_angle_rad << ","
                << n.d_initial << ","
                << n.d_converged << ","
                << n.search_steps << ","
                << n.adjustment << ","
                << n.state << ","
                << n.q.a << "," << n.q.b << "," << n.q.c << "," << n.q.d << "\n";

        std::cout << std::fixed
                  << n.index << ","
                  << std::setprecision(15) << n.riemann_zero_gamma << ","
                  << n.prime << ","
                  << std::setprecision(6) << n.phase_angle_rad << ","
                  << n.d_initial << ","
                  << n.d_converged << ","
                  << n.search_steps << ","
                  << n.adjustment << ","
                  << n.state << ","
                  << n.q.a << "," << n.q.b << "," << n.q.c << "," << n.q.d << "\n";
    }

    return 0;
}
