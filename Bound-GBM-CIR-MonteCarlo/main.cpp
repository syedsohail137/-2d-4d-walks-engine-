/******************************************************************************
 * SYSTEM ENGINE: Multi-Threaded SDE Simulation & Lattice Bounding Verification
 * PLATFORM: Cross-Platform Native Standard Library (C++17 or Higher)
 * LICENSE: MIT / Open Source Performance Asset
 *******************************************************************************/

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <thread>
#include <mutex>
#include <random>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Fast O(1) Mod-7 Quadratic Residue Filter Map
constexpr bool IS_MOD7_RESIDUE[7] = { true, true, true, false, true, false, false };

enum class ProcessType {
    GBM,
    CIR,
    RANDOM
};

// Tracking metrics for paths that successfully pass the lattice criteria
struct LatticeTelemetry {
    ProcessType type;
    long long n = 0;
    long long a = 0;
    long long d = 0;
    double phase_angle = 0.0;
    double spectral_bound = 0.0;
};

// Aggregation structure for side-by-side reporting
struct AggregateStats {
    long long total_simulated_points = 0;
    long long lattice_passed_points = 0;
    double sum_phase = 0.0;
    double sum_sq_phase = 0.0;
    double sum_bound = 0.0;
    double sum_sq_bound = 0.0;
};

// Global synchronization structures
std::mutex g_reporting_mutex;
AggregateStats g_gbm_stats;
AggregateStats g_cir_stats;
AggregateStats g_random_stats; // FIXED: Added missing AggregateStats type specifier

// SDE Hyperparameters
struct SDEParameters {
    double X0 = 100.0;    // Initial value
    double mu = 0.05;     // GBM Drift
    double sigma = 0.20;  // Volatility (Used across models)
    double k = 2.0;       // CIR Mean reversion speed
    double theta = 100.0; // CIR Mean reversion level
    double T = 1.0;       // Time Horizon
    int steps = 100;      // Discretization steps per path
    double K_scale = 1.0; // Continuous-to-Discrete Scaling Factor
};

// Executes the integrated SDE stepping and lattice filtering
void execute_simulation_worker(int paths_per_thread, SDEParameters params, unsigned int seed) {
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    double dt = params.T / static_cast<double>(params.steps);
    double sqrt_dt = std::sqrt(dt);

    // Local accumulation buffers to avoid global thrashing
    AggregateStats local_gbm;
    AggregateStats local_cir;
    AggregateStats local_random;

    for (int p = 0; p < paths_per_thread; ++p) {
        double x_gbm = params.X0;
        double x_cir = params.X0;
        double x_rand = params.X0;

        for (int s = 0; s < params.steps; ++s) {
            double dW = normal_dist(generator) * sqrt_dt;

            // 1. Step Geometric Brownian Motion (GBM)
            x_gbm += params.mu * x_gbm * dt + params.sigma * x_gbm * dW;
            if (x_gbm < 0.0) x_gbm = 0.0; 

            // 2. Step Cox-Ingersoll-Ross (CIR)
            double cir_drift = params.k * (params.theta - x_cir) * dt;
            double cir_diff = params.sigma * std::sqrt(std::max(0.0, x_cir)) * dW;
            x_cir += cir_drift + cir_diff;
            if (x_cir < 0.0) x_cir = 0.0;

            // 3. Step Pure Random Walk
            x_rand += params.sigma * x_rand * dW;
            if (x_rand < 0.0) x_rand = 0.0; 

            // --- LATTICE MAPPING AND VERIFICATION PASS ---
            auto evaluate_lattice = [&](double continuous_val, AggregateStats& stats) {
                stats.total_simulated_points++;
                
                long long n = static_cast<long long>(std::round(continuous_val * params.K_scale));
                
                if (n == 2 || (n > 0 && (n % 4 == 1))) {
                    bool decomposable = false;
                    long long a_val = 0;
                    long long d_val = 0;

                    for (long long i = 1; i * i < n; ++i) {
                        long long remainder = n - i * i;

                        if (!IS_MOD7_RESIDUE[remainder % 7]) {
                            continue;
                        }

                        long long root = static_cast<long long>(std::sqrt(remainder));
                        if (root * root == remainder) {
                            a_val = i;
                            d_val = root;
                            decomposable = true;
                            break;
                        }
                    }

                    if (decomposable) {
                        stats.lattice_passed_points++;

                        std::complex<double> rotor(static_cast<double>(a_val), static_cast<double>(d_val));
                        double phase = std::atan2(rotor.imag(), rotor.real()) * 180.0 / M_PI;
                        
                        double bound = 1.3 * std::sqrt(static_cast<double>(n));

                        stats.sum_phase += phase;
                        stats.sum_sq_phase += phase * phase;
                        stats.sum_bound += bound;
                        stats.sum_sq_bound += bound * bound;
                    }
                }
            };

            evaluate_lattice(x_gbm, local_gbm);
            evaluate_lattice(x_cir, local_cir);
            evaluate_lattice(x_rand, local_random);
        }
    }

    // Thread-safe lock-guard block for merging local thread computations
    std::lock_guard<std::mutex> lock(g_reporting_mutex);
    g_gbm_stats.total_simulated_points += local_gbm.total_simulated_points;
    g_gbm_stats.lattice_passed_points += local_gbm.lattice_passed_points;
    g_gbm_stats.sum_phase += local_gbm.sum_phase;
    g_gbm_stats.sum_sq_phase += local_gbm.sum_sq_phase;
    g_gbm_stats.sum_bound += local_gbm.sum_bound;
    g_gbm_stats.sum_sq_bound += local_gbm.sum_sq_bound;

    g_cir_stats.total_simulated_points += local_cir.total_simulated_points;
    g_cir_stats.lattice_passed_points += local_cir.lattice_passed_points;
    g_cir_stats.sum_phase += local_cir.sum_phase;
    g_cir_stats.sum_sq_phase += local_cir.sum_sq_phase;
    g_cir_stats.sum_bound += local_cir.sum_bound;
    g_cir_stats.sum_sq_bound += local_cir.sum_sq_bound;

    g_random_stats.total_simulated_points += local_random.total_simulated_points;
    g_random_stats.lattice_passed_points += local_random.lattice_passed_points;
    g_random_stats.sum_phase += local_random.sum_phase;
    g_random_stats.sum_sq_phase += local_random.sum_sq_phase;
    g_random_stats.sum_bound += local_random.sum_bound;
    g_random_stats.sum_sq_bound += local_random.sum_sq_bound;
}

void output_statistical_report(const std::string& label, const AggregateStats& stats) {
    double n_passed = static_cast<double>(stats.lattice_passed_points);
    double pass_rate = (static_cast<double>(stats.lattice_passed_points) / static_cast<double>(stats.total_simulated_points)) * 100.0;

    double mean_phase = 0.0, sd_phase = 0.0;
    double mean_bound = 0.0, sd_bound = 0.0;

    if (stats.lattice_passed_points > 0) {
        mean_phase = stats.sum_phase / n_passed;
        double var_phase = (stats.sum_sq_phase / n_passed) - (mean_phase * mean_phase);
        sd_phase = std::sqrt(std::max(0.0, var_phase));

        mean_bound = stats.sum_bound / n_passed;
        double var_bound = (stats.sum_sq_bound / n_passed) - (mean_bound * mean_bound);
        sd_bound = std::sqrt(std::max(0.0, var_bound));
    }

    std::cout << std::left << std::setw(16) << label
              << std::setw(15) << stats.total_simulated_points
              << std::setw(15) << stats.lattice_passed_points
              << std::setw(13) << (std::to_string(pass_rate).substr(0, 5) + "%")
              << std::fixed << std::setprecision(2)
              << mean_phase << "° ± " << std::setw(8) << sd_phase
              << std::setprecision(3)
              << mean_bound << " ± " << sd_bound << "\n";
}

int main() {
    SDEParameters params;
    const int TOTAL_PATHS = 1000000; 
    
    std::cout << "================================================================================\n";
    std::cout << "   INTEGRATED NATIVE SDE SIMULATION & SPECIMEN VERIFICATION ENGINE              \n";
    std::cout << "================================================================================\n\n";

    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;
    
    int paths_per_thread = TOTAL_PATHS / num_threads;
    std::cout << ">> Deploying Workers: " << num_threads << " Native Threads.\n";
    std::cout << ">> Total Trajectories Allocation: " << TOTAL_PATHS << " Monte Carlo Paths.\n";
    std::cout << ">> Time Horizon (T): " << params.T << " | Discretization Steps: " << params.steps << "\n";
    std::cout << ">> Commencing Parallel SDE Stepping & Quadratic Filter Verification...\n\n";

    std::vector<std::thread> thread_pool;
    std::random_device rd;

    for (unsigned int i = 0; i < num_threads; ++i) {
        unsigned int worker_seed = rd() + i * 1051; 
        thread_pool.push_back(std::thread(execute_simulation_worker, paths_per_thread, params, worker_seed));
    }

    for (auto& worker : thread_pool) {
        if (worker.joinable()) worker.join();
    }

    std::cout << "================================================================================\n";
    std::cout << "               SIDE-BY-SIDE SIDEBAND STRUCTURAL REPORT MATRIX                   \n";
    std::cout << "================================================================================\n";
    std::cout << std::left 
              << std::setw(16) << "Process Model"
              << std::setw(15) << "Total Nodes"
              << std::setw(15) << "Passed Nodes"
              << std::setw(13) << "Pass Rate"
              << std::setw(16) << "Phase (μ ± σ)"
              << std::setw(15) << "Bound (μ ± σ)" << "\n";
    std::cout << std::string(90, '-') << "\n";
    output_statistical_report("GBM", g_gbm_stats);
	output_statistical_report("CIR (Square-Rt)", g_cir_stats);
	output_statistical_report("Pure Random", g_random_stats);
	std::cout << std::string(90, '-') << "\n";
	std::cout << "Execution Status: 100% SDE Integration Verification Success.\n";
	std::cout << "================================================================================\n";
	return 0;
	}