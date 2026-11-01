#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdint>
#include <algorithm>
#include <string>
#include <iomanip> // Fixes the std::setprecision / std::fixed error
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

    // Computes the Euclidean squared norm of the accumulated vector state
    double norm_sq() const {
        return static_cast<double>(x * x + y * y);
    }

    // Direct translation vector injection along the bivector plane orientation
    void add_element(int16_t dx, int16_t dy) {
        x += dx;
        y += dy;
    }
};

// High-speed Fermat Decomposer checking sums of two squares
bool fast_decompose_two_square(uint64_t p, int16_t& x, int16_t& y) {
    if (p == 2) { x = 1; y = 1; return true; }
    if ((p & 3) != 1) return false; // Strict filter for primes ≡ 1 (mod 4)

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
    
    // Construct the active Fermat filter subset
    std::vector<uint64_t> test_primes;
    for (auto p : base_primes) {
        if (p == 2 || (p % 4 == 1)) {
            test_primes.push_back(p);
        }
    }

    size_t total_primes = test_primes.size();
    
    // Memory mapping matrix allocations for telemetry tracking tracks
    std::vector<double> res_track_A(total_primes, 0.0);
    std::vector<double> res_track_B(total_primes, 0.0);
    std::vector<double> results_bound(total_primes, 0.0);

    std::vector<int16_t> comp_x(total_primes, 0), comp_y(total_primes, 0);
    std::vector<bool> decomposable(total_primes, false);

    std::cout << "[SYSTEM] Phase 1: Processing parallelized 2D Fermat decompositions..." << std::endl;

    // Distribute coordinate decomposition workload across active CPU cores
    #pragma omp parallel for schedule(static)
    for (size_t idx = 0; idx < total_primes; ++idx) {
        decomposable[idx] = fast_decompose_two_square(test_primes[idx], comp_x[idx], comp_y[idx]);
    }

    std::cout << "[SYSTEM] Phase 2: Running deterministic tracking loop..." << std::endl;

    // Fixed index mapping tables for D4 dihedral transformations
    const int perm_patterns[2][2] = { {0, 1}, {1, 0} };
    const int sign_patterns[4][2] = { {1, 1}, {1, -1}, {-1, 1}, {-1, -1} };

    Clifford2D track_A; // Raw drift tracker
    Clifford2D track_B; // Bivector-stabilized tracker

    for (size_t idx = 0; idx < total_primes; ++idx) {
        if (!decomposable[idx]) continue;

        uint64_t p = test_primes[idx];
        int16_t source_coords[2] = { comp_x[idx], comp_y[idx] };

        // Track A: Baseline raw growth accumulation
        track_A.add_element(source_coords[0], source_coords[1]);

        // Track B: Scan the 8 orbit permutations permitted by the Clifford plane
        double b_best_norm_sq = 1e300;
        int16_t b_chosen_dx = 0, b_chosen_dy = 0;

        for (int p_idx = 0; p_idx < 2; ++p_idx) {
            int64_t cx = source_coords[perm_patterns[p_idx][0]];
            int64_t cy = source_coords[perm_patterns[p_idx][1]];

            for (int s_idx = 0; s_idx < 4; ++s_idx) {
                int16_t sx = static_cast<int16_t>(cx * sign_patterns[s_idx][0]);
                int16_t sy = static_cast<int16_t>(cy * sign_patterns[s_idx][1]);

                Clifford2D candidate_state = track_B;
                candidate_state.add_element(sx, sy);
                double candidate_norm_sq = candidate_state.norm_sq();

                if (candidate_norm_sq < b_best_norm_sq) {
                    b_best_norm_sq = candidate_norm_sq;
                    b_chosen_dx = sx;
                    b_chosen_dy = sy;
                }
            }
        }
        // Commit the stabilizing step that forces destructive interference
        track_B.add_element(b_chosen_dx, b_chosen_dy);

        // Extract and map floating-point matrix trace summaries
        res_track_A[idx] = std::sqrt(track_A.norm_sq());
        res_track_B[idx] = std::sqrt(track_B.norm_sq());
        results_bound[idx] = 1.3 * std::sqrt(static_cast<double>(p));
    }

    // Write the output file using standard telemetry formats
    std::ofstream csv_file("output/fermat_bivector_comparative_res.csv");
    csv_file << "Prime,Track_A_Raw,Track_B_Greedy,Spectral_Bound\n";
    csv_file << std::fixed << std::setprecision(5);
    for (size_t i = 0; i < total_primes; ++i) {
        if (test_primes[i] == 0) continue;
        csv_file << test_primes[i] << "," 
                 << res_track_A[i] << "," 
                 << res_track_B[i] << "," 
                 << results_bound[i] << "\n";
    }
    csv_file.close();

    std::cout << "\n=========================================================" << std::endl;
    std::cout << "[COMPLETE] Metrics compiled to output/fermat_bivector_comparative_res.csv." << std::endl;
    std::cout << "[ANALYSIS] Primes Tracked: " << total_primes << std::endl;
    std::cout << "[ANALYSIS] Final Track A Magnitude: " << res_track_A.back() << std::endl;
    std::cout << "[ANALYSIS] Final Track B Balanced Magnitude: " << res_track_B.back() << std::endl;
    std::cout << "[ANALYSIS] Upper Spectral Boundary Target: " << results_bound.back() << std::endl;
    std::cout << "=========================================================" << std::endl;

    return 0;
}
