#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdint>
#include <algorithm>
#include <string>
#include <iomanip>
#include <limits>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir(dir, mode) _mkdir(dir)
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// High-performance 2D Clifford Element Structure
struct Clifford2D {
    int64_t x = 0;
    int64_t y = 0;

    double norm_sq() const {
        return static_cast<double>(x * x + y * y);
    }
};

// Computes the Logarithmic Barrier Penalty to enforce Wolfovitz-Vector Balancing
inline double compute_log_barrier_2d(int64_t x, int64_t y, double boundary) {
    double b_sq = boundary * boundary;
    double total_penalty = 0.0;
    
    double coords[2] = { static_cast<double>(x), static_cast<double>(y) };
    
    for (int i = 0; i < 2; ++i) {
        double margin = b_sq - (coords[i] * coords[i]);
        if (margin <= 1e-5) {
            return std::numeric_limits<double>::infinity(); // Hard wall barrier ceiling limit reached
        }
        total_penalty -= std::log(margin); // Symmetric Log-Barrier cost evaluation
    }
    return total_penalty;
}

// High-speed Fermat Decomposer checking sums of two squares
bool fast_decompose_two_square(uint64_t p, int16_t& x, int16_t& y) {
    if (p == 2) { x = 1; y = 1; return true; }
    if ((p & 3) != 1) return false; 

    uint64_t limit_x = static_cast<uint64_t>(std::sqrt(p));
    for (uint64_t tx = limit_x; tx > 0; --tx) {
        uint64_t rem = p - (tx * tx);
        uint64_t ty = static_cast<uint64_t>(std::sqrt(rem));
        if (tx * tx + ty * ty == p) {
            x = static_cast<int16_t>(tx);
            y = static_cast<int16_t>(ty);
            return true;
        }
    }
    return false;
}

// Optimized Linear Sieve for generation of base prime matrices
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
    // Generate output structure paths safely on Posix or Windows
    mkdir("output", 0777);

    uint64_t upper_bound_limit = 100000; 
    std::vector<uint64_t> base_primes = generate_high_range_primes(upper_bound_limit);
    
    std::vector<uint64_t> test_primes;
    for (auto p : base_primes) {
        if (p == 2 || (p % 4 == 1)) {
            test_primes.push_back(p);
        }
    }

    size_t total_primes = test_primes.size();
    
    std::vector<double> res_track_A(total_primes, 0.0);
    std::vector<double> res_track_B(total_primes, 0.0);
    std::vector<double> results_bound(total_primes, 0.0);

    std::vector<int16_t> comp_x(total_primes, 0), comp_y(total_primes, 0);
    std::vector<bool> decomposable(total_primes, false);

    std::cout << "[SYSTEM] Step 1: Processing parallelized 2D Fermat decompositions..." << std::endl;

    #pragma omp parallel for schedule(static)
    for (size_t idx = 0; idx < total_primes; ++idx) {
        decomposable[idx] = fast_decompose_two_square(test_primes[idx], comp_x[idx], comp_y[idx]);
    }

    std::cout << "[SYSTEM] Step 2: Running non-exponential log-barrier tracking loop..." << std::endl;

    const int perm_patterns[2][2] = { {0, 1}, {1, 0} };

    Clifford2D track_A; // Raw baseline tracking
    Clifford2D track_B; // Structural logarithmic barrier tracking

    for (size_t idx = 0; idx < total_primes; ++idx) {
        if (!decomposable[idx]) continue;

        uint64_t p = test_primes[idx];
        double current_bound = 1.3 * std::sqrt(static_cast<double>(p));
        int16_t source_coords[2] = { comp_x[idx], comp_y[idx] };

        // Track A: Baseline raw growth accumulation
        track_A.x += source_coords[0];
        track_A.y += source_coords[1];

        // Track B: Scan permutations evaluating localized log barrier fields
        double best_barrier_score = std::numeric_limits<double>::infinity();
        int16_t b_chosen_dx = 0, b_chosen_dy = 0;

        for (int p_idx = 0; p_idx < 2; ++p_idx) {
            int64_t cx = source_coords[perm_patterns[p_idx][0]];
            int64_t cy = source_coords[perm_patterns[p_idx][1]];

            // O(1) Analytic Sign Resolution: Forces direct destructive counter-drift
            int16_t sx = static_cast<int16_t>((track_B.x >= 0) ? -cx : cx);
            int16_t sy = static_cast<int16_t>((track_B.y >= 0) ? -cy : cy);

            double score = compute_log_barrier_2d(track_B.x + sx, track_B.y + sy, current_bound);

            if (score < best_barrier_score) {
                best_barrier_score = score;
                b_chosen_dx = sx;
                b_chosen_dy = sy;
            }
        }

        // Fallback Scheme: If trapped inside an asymptotic singularity, fall back to pure local greedy steps
        if (best_barrier_score == std::numeric_limits<double>::infinity()) {
            int16_t cx = source_coords[0]; int16_t cy = source_coords[1];
            b_chosen_dx = static_cast<int16_t>((track_B.x >= 0) ? -cx : cx);
            b_chosen_dy = static_cast<int16_t>((track_B.y >= 0) ? -cy : cy);
        }

        // Commit the stabilizing step
        track_B.x += b_chosen_dx;
        track_B.y += b_chosen_dy;

        // Map magnitudes to arrays
        res_track_A[idx] = std::sqrt(track_A.norm_sq());
        res_track_B[idx] = std::sqrt(track_B.norm_sq());
        results_bound[idx] = current_bound;
    }

    // Write out metrics to comparative matrix logs
    std::ofstream csv_file("output/fermat_barrier_comparative_res.csv");
    csv_file << "Prime,Track_A_Raw,Track_B_Log_Barrier,Spectral_Bound\n";
    csv_file << std::fixed << std::setprecision(5);
    for (size_t i = 0; i < total_primes; ++i) {
        if (test_primes[i] == 0) continue;
        csv_file << test_primes[i] << "," 
                 << res_track_A[i] << "," 
                 << res_track_B[i] << "," 
                 << results_bound[i] << "\n";
    }
    csv_file.close();

    std::cout << "[COMPLETE] Metrics compiled to: output/fermat_barrier_comparative_res.csv" << std::endl;
    return 0;
}
