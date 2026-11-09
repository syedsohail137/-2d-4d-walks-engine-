#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>

// Re-packing data structs for multi-phasor tracking
struct ExtendedPrimeNode4D {
    uint64_t prime;
    int16_t a, b, c, d;
    double cumulative_phase_r;
    double cumulative_phase_i;
};

struct MertensDrift4D {
    double cumulative_x = 0.0;
    double cumulative_y = 0.0;
    double cumulative_z = 0.0;
    double cumulative_w = 0.0;

    double compute_net_displacement() const {
        return std::sqrt(cumulative_x * cumulative_x + 
                         cumulative_y * cumulative_y + 
                         cumulative_z * cumulative_z + 
                         cumulative_w * cumulative_w);
    }
};

// Precise 4-Square decomposition matching your blueprint constraints
bool decompose_four_square(uint64_t p, int16_t& a, int16_t& b, int16_t& c, int16_t& d) {
    int64_t max_a = static_cast<int64_t>(std::sqrt(p));
    for (int64_t ta = max_a; ta >= 0; --ta) {
        int64_t rem_a = p - (ta * ta); if (rem_a < 0) continue;
        int64_t max_b = static_cast<int64_t>(std::sqrt(rem_a));
        for (int64_t tb = max_b; tb >= 0; --tb) {
            int64_t rem_b = rem_a - (tb * tb); if (rem_b < 0) continue;
            int64_t max_c = static_cast<int64_t>(std::sqrt(rem_b));
            for (int64_t tc = max_c; tc >= 0; --tc) {
                int64_t rem_c = rem_b - (tc * tc);
                int64_t td = static_cast<int64_t>(std::round(std::sqrt(rem_c)));
                
                if (ta * ta + tb * tb + tc * tc + td * td == p) {
                    a = static_cast<int16_t>(ta); b = static_cast<int16_t>(tb);
                    c = static_cast<int16_t>(tc); d = static_cast<int16_t>(td);
                    return true;
                }
            }
        }
    }
    return false;
}

int main() {
    // Structural frequencies configuration array
    const std::vector<double> zeta_gammas = {
        14.134725141734693,
        21.022039638771555,
        25.010857580145688
    };

    std::vector<uint64_t> test_primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
    MertensDrift4D running_drift;
    
    std::cout << "[MULTI-PHASOR SIMULATION] Initiating Interference Check (γ₁ + γ₂ + γ₃)...\n";
    std::cout << "------------------------------------------------------------------------\n";
    std::cout << std::setw(6) << "Prime" << " | " 
              << std::setw(15) << "Coordinates" << " | " 
              << std::setw(15) << "Net Drift" << " | " 
              << std::setw(15) << "Riemann Bound" << " | Status\n";
    std::cout << "------------------------------------------------------------------------\n";

    uint64_t current_N = 0;
    for (size_t i = 0; i < test_primes.size(); ++i) {
        uint64_t p = test_primes[i];
        current_N = p; // Evaluating bound at the current upper prime limit
        
        int16_t a = 0, b = 0, c = 0, d = 0;
        if (decompose_four_square(p, a, b, c, d)) {
            double log_p = std::log(static_cast<double>(p));
            double sum_cos = 0.0, sum_sin = 0.0;

            // Combine phasors linearly to simulate structural acoustic vibrations
            for (double gamma : zeta_gammas) {
                sum_cos += std::cos(gamma * log_p);
                sum_sin += std::sin(gamma * log_p);
            }

            // Average the superimposed phasor fields
            double phase_r = sum_cos / zeta_gammas.size();
            double phase_i = sum_sin / zeta_gammas.size();

            // Accumulate global spatial walk
            running_drift.cumulative_x += a * phase_r;
            running_drift.cumulative_y += b * phase_i;
            running_drift.cumulative_z += c * phase_r;
            running_drift.cumulative_w += d * phase_i;

            double net_drift = running_drift.compute_net_displacement();
            double riemann_bound = std::pow(static_cast<double>(current_N), 0.5); // O(N^0.5) baseline

            std::string coordinates_str = "(" + std::to_string(a) + "," + std::to_string(b) + "," + std::to_string(c) + "," + std::to_string(d) + ")";
            
            std::cout << std::setw(6) << p << " | " 
                      << std::setw(15) << coordinates_str << " | " 
                      << std::setw(15) << std::fixed << std::setprecision(4) << net_drift << " | " 
                      << std::setw(15) << riemann_bound << " | "
                      << (net_drift <= riemann_bound ? "SAFE" : "OVER LIMIT") << "\n";
        }
    }
    return 0;
}
