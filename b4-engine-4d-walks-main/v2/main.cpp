#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <string>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

// 4D COMPACT METRIC TRACKER (CLOSED TOPOLOGY MOTOR)
struct Compact4DPhaseCage {
    int64_t phase_accumulator_r = 0;
    int64_t phase_accumulator_i = 0;
    const int64_t phase_scale_base = 10000000;

    void step_curved_phasor(uint64_t prime, const std::vector<double>& gammas, int64_t braid_twist) {
        double log_p = std::log(static_cast<double>(prime));
        double sum_cos = 0.0;
        double sum_sin = 0.0;

        // Dynamic spatial throttling: transforms torque into internal metric curvature
        double dynamic_damping = 1.0 / (1.0 + std::log(1.0 + std::abs(braid_twist)));
        double anticommutative_skew = static_cast<double>(braid_twist) * 0.00002 * dynamic_damping;

        #pragma omp simd reduction(+:sum_cos, sum_sin)
        for (size_t k = 0; k < gammas.size(); ++k) {
            double localized_angle = (gammas[k] * log_p) + anticommutative_skew;
            sum_cos += std::cos(localized_angle);
            sum_sin += std::sin(localized_angle);
        }

        phase_accumulator_r = static_cast<int64_t>((sum_cos / gammas.size()) * phase_scale_base);
        phase_accumulator_i = static_cast<int64_t>((sum_sin / gammas.size()) * phase_scale_base);
    }
};

struct BraidCrossing4D {
    int64_t twist_index = 0;         
    int64_t crossing_over_under = 0; 
    
    void process_factor_braid(int16_t a, int16_t b, int16_t c, int16_t d) {
        int64_t orientation_determinant = (a * b) - (c * d);
        if (orientation_determinant % 2 == 0) {
            crossing_over_under += 1;
        } else {
            crossing_over_under -= 1;
        }
        twist_index += (a * d) - (b * c);
    }
};

inline bool decompose_four_square(uint64_t p, int16_t& a, int16_t& b, int16_t& c, int16_t& d) {
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

std::vector<uint64_t> generate_primes_to_limit(uint64_t limit) {
    std::vector<bool> is_prime(limit + 1, true);
    std::vector<uint64_t> primes;
    if (limit < 2) return primes;
    for (uint64_t p = 2; p * p <= limit; ++p) {
        if (is_prime[p]) {
            for (uint64_t i = p * p; i <= limit; i += p) is_prime[i] = false;
        }
    }
    for (uint64_t p = 2; p <= limit; ++p) {
        if (is_prime[p]) primes.push_back(p);
    }
    return primes;
}

void run_benchmark(uint64_t horizon_limit, const std::vector<double>& gammas) {
    std::vector<uint64_t> test_primes = generate_primes_to_limit(horizon_limit);
    size_t total_primes = test_primes.size();

    double total_quantum_x = 0.0, total_quantum_y = 0.0, total_quantum_z = 0.0, total_quantum_w = 0.0;
    int64_t global_twist_index = 0;
    uint64_t critical_shifts_encountered = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for reduction(+:total_quantum_x, total_quantum_y, total_quantum_z, total_quantum_w, global_twist_index, critical_shifts_encountered) schedule(dynamic, 64)
    for (size_t i = 0; i < total_primes; ++i) {
        uint64_t p = test_primes[i];
        int16_t a = 0, b = 0, c = 0, d = 0;

        if (decompose_four_square(p, a, b, c, d)) {
            Compact4DPhaseCage local_cage;
            BraidCrossing4D local_braid;

            local_braid.process_factor_braid(a, b, c, d);
            local_cage.step_curved_phasor(p, gammas, local_braid.twist_index);

            double normalized_r = static_cast<double>(local_cage.phase_accumulator_r) / local_cage.phase_scale_base;
            double normalized_i = static_cast<double>(local_cage.phase_accumulator_i) / local_cage.phase_scale_base;
            
            double horizon_scale_damping = 1.0 / std::log(static_cast<double>(p));
            double topological_correction = (1.0 + (static_cast<double>(local_braid.crossing_over_under) * 0.05)) * horizon_scale_damping;

            double step_x = (a * normalized_r) * topological_correction;
            double step_y = (b * normalized_i) * topological_correction;
            double step_z = (c * normalized_r) * topological_correction;
            double step_w = (d * normalized_i) * topological_correction;

            total_quantum_x += step_x;
            total_quantum_y += step_y;
            total_quantum_z += step_z;
            total_quantum_w += step_w;
            global_twist_index += local_braid.twist_index;

            double local_drift = std::sqrt(step_x*step_x + step_y*step_y + step_z*step_z + step_w*step_w);
            if (local_drift > std::pow(static_cast<double>(p), 0.5)) {
                critical_shifts_encountered++;
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> execution_ms = end_time - start_time;

    double final_net_displacement = std::sqrt(total_quantum_x * total_quantum_x + total_quantum_y * total_quantum_y + 
                                             total_quantum_z * total_quantum_z + total_quantum_w * total_quantum_w);

    std::cout << " > Horizon Limit (N)       : " << horizon_limit << "\n";
    std::cout << " > Execution Runtime        : " << execution_ms.count() << " ms\n";
    std::cout << " > Global Net Trajectory Drift: " << final_net_displacement << "\n";
    std::cout << " > Ultimate Riemann Bound   : " << std::sqrt(horizon_limit) << "\n";
    std::cout << " > Global Braid Twist Index : " << global_twist_index << "\n";
    std::cout << " > Critical Shift Gaps      : " << critical_shifts_encountered << " nodes leaked through\n";
    std::cout << "----------------------------------------------------------------------------------------\n";
}

int main() {
    const std::vector<double> resonant_wall_gammas = {
        14.1347251417, 21.0220396388, 37.5861781588, 52.9703214777,
        75.7046906991, 92.4918992705, 111.8746591771, 143.1118458076
    };

    std::cout << "[INIT] Running Benchmarks for Scaling Horizons...\n";
    std::cout << "----------------------------------------------------------------------------------------\n";
    run_benchmark(100000, resonant_wall_gammas);
    run_benchmark(1000000, resonant_wall_gammas);
    
    // 100M horizon placeholder run simulation note
    std::cout << "[NOTE] For N=100M, run on an enterprise multicore server environment.\n";
    return 0;
}
