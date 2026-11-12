#ifndef IO_AND_MONTECARLO_HPP
#define IO_AND_MONTECARLO_HPP

#include "lattice_solver.hpp"
#include "spectral_inversion.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace spectral {

// Parses external LMFDB or Odlyzko zero tables (supports "index gamma" or single-column formats)
inline std::vector<double> load_gammas_from_file(const std::string& filepath, std::size_t max_count) {
    std::ifstream in(filepath);
    if (!in) {
        throw std::runtime_error("Could not open gammas file: " + filepath);
    }

    std::vector<double> gammas;
    std::string line;
    while (std::getline(in, line)) {
        std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#' || line[first] == '/') {
            continue;
        }

        std::istringstream iss(line.substr(first));
        std::vector<double> tokens;
        double val = 0.0;
        while (iss >> val) {
            tokens.push_back(val);
        }

        if (tokens.empty()) continue;

        double gamma = (tokens.size() >= 2 && tokens[0] == std::floor(tokens[0]) && tokens[1] > 10.0)
		? tokens[1]
	: tokens[0];


        if (gamma > 0.0) {
            gammas.push_back(gamma);
            if (max_count > 0 && gammas.size() == max_count) break;
        }
    }

    if (gammas.empty()) {
        throw std::runtime_error("No valid positive Riemann zero values parsed from: " + filepath);
    }
    return gammas;
}

// Runs multi-seed Monte Carlo trials using uniform phase offsets in [0, 2*PI)
inline MonteCarloSummary run_monte_carlo_comparison(
    const std::vector<SpectralNode>& riemann_nodes,
    long long max_steps,
    int num_trials,
    uint64_t seed) {

    MonteCarloSummary stats;
    stats.trials = num_trials;
    stats.node_count = riemann_nodes.size();

    long long riemann_steps_sum = 0;
    std::size_t riemann_solved = 0;
    for (const auto& n : riemann_nodes) {
        if (n.state == "Exact Match") ++stats.riemann_exact_matches;
        if (n.state != "No Representation") {
            riemann_steps_sum += n.search_steps;
            ++riemann_solved;
        }
    }
    stats.riemann_mean_steps = riemann_solved > 0
        ? static_cast<double>(riemann_steps_sum) / riemann_solved
        : 0.0;

    if (num_trials <= 0 || riemann_nodes.empty()) return stats;

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> phase_dist(0.0, 2.0 * PI);

    double sum_exact = 0.0;
    double sum_exact_sq = 0.0;
    double sum_mean_steps = 0.0;
    int ge_riemann_count = 0;

    for (int t = 0; t < num_trials; ++t) {
        std::size_t trial_exact = 0;
        long long trial_steps = 0;
        std::size_t trial_solved = 0;

        for (std::size_t i = 0; i < riemann_nodes.size(); ++i) {
            long long p = riemann_nodes[i].prime;
            double rand_theta = phase_dist(rng);
            long long rand_d_init = std::llround(static_cast<double>(p) * 0.5 * std::cos(rand_theta));

            SpectralNode sim_node = resolve_from_initial_d(
                i, riemann_nodes[i].riemann_zero_gamma, p, rand_theta, rand_d_init, max_steps);

            if (sim_node.state == "Exact Match") ++trial_exact;
            if (sim_node.state != "No Representation") {
                trial_steps += sim_node.search_steps;
                ++trial_solved;
            }
        }

        double exact_d = static_cast<double>(trial_exact);
        sum_exact += exact_d;
        sum_exact_sq += exact_d * exact_d;
        if (trial_solved > 0) {
            sum_mean_steps += static_cast<double>(trial_steps) / trial_solved;
        }
        if (trial_exact >= stats.riemann_exact_matches) {
            ++ge_riemann_count;
        }
    }

    stats.baseline_mean_exact = sum_exact / num_trials;
    double variance = (sum_exact_sq / num_trials) - (stats.baseline_mean_exact * stats.baseline_mean_exact);
    stats.baseline_std_exact = variance > 0.0 ? std::sqrt(variance) : 0.0;
    stats.baseline_mean_steps = sum_mean_steps / num_trials;

    stats.z_score_exact = (stats.baseline_std_exact > 1e-12)
        ? (static_cast<double>(stats.riemann_exact_matches) - stats.baseline_mean_exact) / stats.baseline_std_exact
        : 0.0;

    stats.empirical_p_value = static_cast<double>(ge_riemann_count + 1) / (num_trials + 1);
    return stats;
}

// Writes combined side-by-side CSV with both Z^4 lattice coordinates and exact Q rational lifts
inline bool export_side_by_side_csv(
    const std::string& filepath,
    const std::vector<SpectralNode>& riemann_nodes,
    const std::vector<SpectralNode>& baseline_nodes) {

    std::ofstream out(filepath);
    if (!out) return false;

    out << "Index,Gamma_Full,Prime,"
        << "Riemann_Phase,Riemann_D_Init,Riemann_D_Conv,Riemann_Steps,Riemann_State,"
        << "R_a,R_b,R_c,R_d,R_Q_Det,R_Q_NormSq,R_Q_StereoX,R_Q_StereoY,R_Q_StereoZ,"
        << "Base_Phase,Base_D_Init,Base_D_Conv,Base_Steps,Base_State,"
        << "B_a,B_b,B_c,B_d,B_Q_Det,B_Q_StereoX,B_Q_StereoY,B_Q_StereoZ,"
        << "Step_Diff_Base_Minus_Riemann\n";

    for (std::size_t i = 0; i < riemann_nodes.size(); ++i) {
        const auto& r = riemann_nodes[i];
        const auto& b = baseline_nodes[i];
        int step_diff = b.search_steps - r.search_steps;

        out << std::fixed
            << r.index << ","
            << std::setprecision(15) << r.riemann_zero_gamma << ","
            << r.prime << ","
            << std::setprecision(6) << r.phase_angle_rad << ","
            << r.d_initial << "," << r.d_converged << "," << r.search_steps << "," << r.state << ","
            << r.q.a << "," << r.q.b << "," << r.q.c << "," << r.q.d << ","
            << r.rat.norm_det << "," << r.rat.norm_sq << ","
            << r.rat.stereo_x << "," << r.rat.stereo_y << "," << r.rat.stereo_z << ","
            << std::setprecision(6) << b.phase_angle_rad << ","
            << b.d_initial << "," << b.d_converged << "," << b.search_steps << "," << b.state << ","
            << b.q.a << "," << b.q.b << "," << b.q.c << "," << b.q.d << ","
            << b.rat.norm_det << ","
            << b.rat.stereo_x << "," << b.rat.stereo_y << "," << b.rat.stereo_z << ","
            << step_diff << "\n";
    }
    return true;
}

} // namespace spectral

#endif // IO_AND_MONTECARLO_HPP
