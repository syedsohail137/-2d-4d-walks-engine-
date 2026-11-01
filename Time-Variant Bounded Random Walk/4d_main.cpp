#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdint>
#include <algorithm>
#include <iomanip>
#include <string>
#include <limits>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "head.h"

// High-Performance 4D Vector Space Tracker
struct Clifford4D {
    int64_t e1 = 0; int64_t e2 = 0; int64_t e3 = 0; int64_t e4 = 0;

    double norm_sq_projected(double phase_r, double phase_i) const {
        double d1 = static_cast<double>(e1) * phase_r;
        double d2 = static_cast<double>(e2) * phase_i;
        double d3 = static_cast<double>(e3) * phase_r;
        double d4 = static_cast<double>(e4) * phase_i;
        return (d1*d1 + d2*d2 + d3*d3 + d4*d4);
    }
};

// Computes the Logarithmic Barrier Penalty to enforce Wolfovitz-Vector Balancing
inline double compute_log_barrier(int64_t e1, int64_t e2, int64_t e3, int64_t e4, double boundary) {
    double b_sq = boundary * boundary;
    double coords[4] = {
        static_cast<double>(e1),
        static_cast<double>(e2),
        static_cast<double>(e3),
        static_cast<double>(e4)
    };
    
    double total_penalty = 0.0;
    for (int i = 0; i < 4; ++i) {
        double margin = b_sq - (coords[i] * coords[i]);
        if (margin <= 1e-5) {
            return std::numeric_limits<double>::infinity(); // Immediate hard barrier block
        }
        total_penalty -= std::log(margin); // Symmetric Log-Barrier expansion
    }
    return total_penalty;
}

// Fast Four-Square Decomposer utilizing Legendre Filters
bool fast_decompose_four_square(uint64_t p, int16_t& a, int16_t& b, int16_t& c, int16_t& d) {
    if (p == 2) { a = 1; b = 1; c = 0; d = 0; return true; }
    uint64_t limit_a = static_cast<uint64_t>(std::sqrt(p));
    for (uint64_t ta = limit_a; ta > 0; --ta) {
        uint64_t rem_a = p - (ta * ta);
        if ((rem_a & 7) == 7) continue; 
        uint64_t limit_b = static_cast<uint64_t>(std::sqrt(rem_a));
        for (uint64_t tb = limit_b; tb >= ta; --tb) { 
            uint64_t rem_b = rem_a - (tb * tb);
            uint64_t tc = static_cast<uint64_t>(std::sqrt(rem_b));
            if (tc > 0) {
                uint64_t rem_c = rem_b - (tc * tc);
                uint64_t td = static_cast<uint64_t>(std::sqrt(rem_c));
                if (ta * ta + tb * tb + tc * tc + td * td == p) {
                    a = static_cast<int16_t>(ta); 
                    b = static_cast<int16_t>(tb); 
                    c = static_cast<int16_t>(tc); 
                    d = static_cast<int16_t>(td);
                    return true;
                }
            }
        }
    }
    return false;
}

std::vector<uint64_t> generate_high_range_primes(uint64_t limit) {
    std::vector<bool> is_prime(limit + 1, true);
    std::vector<uint64_t> primes;
    if (limit >= 2) primes.push_back(2);
    for (uint64_t i = 3; i * i <= limit; i += 2) {
        if (is_prime[i]) {
            for (uint64_t j = i * i; j <= limit; j += i * 2) is_prime[j] = false;
        }
    }
    for (uint64_t i = 3; i <= limit; i += 2) {
        if (is_prime[i]) primes.push_back(i);
    }
    return primes;
}

int main() {
#if defined(_WIN32)
    _mkdir("output");
#else
    mkdir("output", 0777);
#endif

    const std::vector<double>& zeta_gammas = SpectralData::zeta_gammas;
    uint64_t upper_bound_limit = 100000; 
    std::vector<uint64_t> test_primes = generate_high_range_primes(upper_bound_limit); 
    size_t total_primes = test_primes.size();
    
    std::vector<double> res_track_A(total_primes, 0.0);
    std::vector<double> res_track_B(total_primes, 0.0);
    std::vector<double> results_bound(total_primes, 0.0);

    std::vector<int16_t> comp_a(total_primes, 0), comp_b(total_primes, 0), comp_c(total_primes, 0), comp_d(total_primes, 0);
    std::vector<double> phase_r_arr(total_primes, 0.0), phase_i_arr(total_primes, 0.0);
    std::vector<bool> decomposable(total_primes, false);

    std::cout << "[SYSTEM] Executing Step 1: Pre-computing Riemann Zeta Transformations..." << std::endl;

    #pragma omp parallel for schedule(static)
    for (size_t idx = 0; idx < total_primes; ++idx) {
        uint64_t p = test_primes[idx];
        double log_p = std::log(static_cast<double>(p));
        double sum_cos = 0.0, sum_sin = 0.0;
        for (double gamma : zeta_gammas) {
            sum_cos += std::cos(gamma * log_p);
            sum_sin += std::sin(gamma * log_p);
        }
        phase_r_arr[idx] = sum_cos / zeta_gammas.size();
        phase_i_arr[idx] = sum_sin / zeta_gammas.size();
        decomposable[idx] = fast_decompose_four_square(p, comp_a[idx], comp_b[idx], comp_c[idx], comp_d[idx]);
    }

    std::cout << "[SYSTEM] Executing Step 2: Non-Exponential Boundary Tracker Engine..." << std::endl;

    // Symmetric 24 Permutations mapping the Weyl orbit lines
    const int perm_patterns[24][4] = {
        {0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1}, {0,3,1,2}, {0,3,2,1},
        {1,0,2,3}, {1,0,3,2}, {1,2,0,3}, {1,2,3,0}, {1,3,0,2}, {1,3,2,0},
        {2,0,1,3}, {2,0,3,1}, {2,1,0,3}, {2,1,3,0}, {2,3,0,1}, {2,3,1,0},
        {3,0,1,2}, {3,0,2,1}, {3,1,0,2}, {3,1,2,0}, {3,2,0,1}, {3,2,1,0}
    };

    Clifford4D track_A; // Raw Control Loop
    Clifford4D track_B; // Logarithmic Barrier Shielded Vector Tracker

    for (size_t idx = 0; idx < total_primes; ++idx) {
        uint64_t p = test_primes[idx];
        double pr = phase_r_arr[idx];
        double pi = phase_i_arr[idx];
        double current_bound = 1.3 * std::sqrt(static_cast<double>(p));

        if (decomposable[idx]) {
            // Track A Tracking Loop (Raw Drift)
            track_A.e1 += comp_a[idx]; track_A.e2 += comp_b[idx]; 
            track_A.e3 += comp_c[idx]; track_A.e4 += comp_d[idx];

            // Track B Logarithmic Barrier Selection Scan
            double best_barrier_score = std::numeric_limits<double>::infinity();
            int16_t b_ea = 0, b_eb = 0, b_ec = 0, b_ed = 0;
            int16_t base_components[4] = { comp_a[idx], comp_b[idx], comp_c[idx], comp_d[idx] };

            for (int p_idx = 0; p_idx < 24; ++p_idx) {
                int64_t cx = base_components[perm_patterns[p_idx][0]];
                int64_t cy = base_components[perm_patterns[p_idx][1]];
                int64_t cz = base_components[perm_patterns[p_idx][2]];
                int64_t cw = base_components[perm_patterns[p_idx][3]];

                // Inject O(1) Analytic Sign Flips relative to structural tracking arrays
                int16_t sx = static_cast<int16_t>((track_B.e1 >= 0) ? -cx : cx);
                int16_t sy = static_cast<int16_t>((track_B.e2 >= 0) ? -cy : cy);
                int16_t sz = static_cast<int16_t>((track_B.e3 >= 0) ? -cz : cz);
                int16_t sw = static_cast<int16_t>((track_B.e4 >= 0) ? -cw : cw);

                // Score candidate step under localized Logarithmic Potential Barrier Field
                double score = compute_log_barrier(track_B.e1 + sx, track_B.e2 + sy, track_B.e3 + sz, track_B.e4 + sw, current_bound);

                if (score < best_barrier_score) {
                    best_barrier_score = score;
                    b_ea = sx; b_eb = sy; b_ec = sz; b_ed = sw;
                }
            }

            // Fallback Recovery Scheme: If a geometric trap is triggered, force raw minimal greedy step
            if (best_barrier_score == std::numeric_limits<double>::infinity()) {
                int16_t cx = comp_a[idx]; int16_t cy = comp_b[idx];
                int16_t cz = comp_c[idx]; int16_t cw = comp_d[idx];
                b_ea = static_cast<int16_t>((track_B.e1 >= 0) ? -cx : cx);
                b_eb = static_cast<int16_t>((track_B.e2 >= 0) ? -cy : cy);
                b_ec = static_cast<int16_t>((track_B.e3 >= 0) ? -cz : cz);
                b_ed = static_cast<int16_t>((track_B.e4 >= 0) ? -cw : cw);
            }

            // Commit the balanced vector update
            track_B.e1 += b_ea; track_B.e2 += b_eb; track_B.e3 += b_ec; track_B.e4 += b_ed;
        }

        res_track_A[idx] = std::sqrt(track_A.norm_sq_projected(pr, pi));
        res_track_B[idx] = std::sqrt(track_B.norm_sq_projected(pr, pi));
        results_bound[idx] = current_bound;
    }

    // Export structural matrix lines to telemetry file
    std::ofstream csv_file("output/wolfovitz_barrier_res.csv");
    csv_file << "Prime,Track_A_Raw,Track_B_Wolfovitz_Barrier,Spectral_Bound\n";
    csv_file << std::fixed << std::setprecision(5);
    for (size_t i = 0; i < total_primes; ++i) {
        csv_file << test_primes[i] << "," << res_track_A[i] << "," << res_track_B[i] << "," << results_bound[i] << "\n";
    }
    csv_file.close();

    std::cout << "[COMPLETE] Matrix logs compiled to: output/wolfovitz_barrier_res.csv" << std::endl;
    return 0;
}
