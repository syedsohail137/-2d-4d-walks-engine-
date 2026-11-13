#include "gue_unfolding.hpp"
#include "io_and_montecarlo.hpp"
#include "lattice_solver.hpp"
#include "phase_accumulator.hpp"
#include "rational_projection.hpp"
#include "spectral_inversion.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "Options:\n"
              << "  --gammas-file <path>  Load Riemann zeros from external LMFDB/Odlyzko file\n"
              << "  --count <N>           Number of prime-zero pairs to evaluate (default: 100)\n"
              << "  --max-steps <N>       Maximum determinant search radius (default: 40)\n"
              << "  --mc-trials <N>       Monte Carlo random phase trials for p-value (default: 500)\n"
              << "  --seed <N>            64-bit RNG seed for baseline & Monte Carlo (default: 42)\n"
              << "  --output <path>       Side-by-side CSV path (default: spectral_comparison_100.csv)\n"
              << "  --verify-gue          Run Module 3 GUE semicircle CDF unfolding check\n"
              << "  --gue-size <N>        Hermitian matrix dimension for GUE check (default: 200)\n"
              << "  --quiet               Suppress per-row stdout preview\n"
              << "  --help                Display this help message\n";
}

bool parse_cli(int argc, char* argv[], spectral::RunConfig& cfg) {
    const char* prog_name = (argc > 0 && argv[0] != nullptr) ? argv[0] : "spectral_inversion";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(prog_name);
            std::exit(0);
        } else if (arg == "--quiet") {
            cfg.quiet = true;
        } else if (arg == "--verify-gue") {
            cfg.verify_gue = true;
        } else if (arg == "--gammas-file" && i + 1 < argc) {
            cfg.gammas_file = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            cfg.output_csv = argv[++i];
        } else if (arg == "--count" && i + 1 < argc) {
            long long val = std::stoll(argv[++i]);
            if (val <= 0) {
                std::cerr << "Error: --count must be positive.\n";
                return false;
            }
            cfg.count = static_cast<std::size_t>(val);
        } else if (arg == "--max-steps" && i + 1 < argc) {
            cfg.max_steps = std::stoll(argv[++i]);
            if (cfg.max_steps < 0) {
                std::cerr << "Error: --max-steps must be non-negative.\n";
                return false;
            }
        } else if (arg == "--mc-trials" && i + 1 < argc) {
            cfg.mc_trials = std::stoi(argv[++i]);
            if (cfg.mc_trials < 0) {
                std::cerr << "Error: --mc-trials must be non-negative.\n";
                return false;
            }
        } else if (arg == "--seed" && i + 1 < argc) {
            cfg.seed = std::stoull(argv[++i]);
        } else if (arg == "--gue-size" && i + 1 < argc) {
            cfg.gue_size = std::stoi(argv[++i]);
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            print_usage(prog_name);
            return false;
        }
    }
    return true;
}

bool verify_node_invariants(const spectral::SpectralNode& n, const char* label) {
    if (n.state != "No Representation") {
        long long norm_sq = n.q.a * n.q.a + n.q.b * n.q.b + n.q.c * n.q.c + n.q.d * n.q.d;
        long long det = n.q.a * n.q.d - n.q.b * n.q.c;
        if (norm_sq != n.prime || det != n.d_converged) {
            std::cerr << "Z^4 Algebraic invariant failure [" << label << "] at Index " << n.index << "\n";
            return false;
        }
        if (!spectral::verify_rational_invariants(n)) {
            std::cerr << "Q Rational field invariant failure [" << label << "] at Index " << n.index << "\n";
            return false;
        }
    }
    return true;
}

} // namespace

int main(int argc, char* argv[]) {
    spectral::RunConfig cfg;
    if (!parse_cli(argc, argv, cfg)) return 1;

    std::vector<double> gammas;
    try {
        if (!cfg.gammas_file.empty()) {
            gammas = spectral::load_gammas_from_file(cfg.gammas_file, cfg.count);
        } else {
            gammas = spectral::DEFAULT_RIEMANN_GAMMAS;
            if (cfg.count > 0 && cfg.count < gammas.size()) {
                gammas.resize(cfg.count);
            } else if (cfg.count > gammas.size()) {
                std::cerr << "Error: Requested --count " << cfg.count
                          << " exceeds embedded table (" << gammas.size()
                          << "). Supply --gammas-file for larger runs.\n";
                return 1;
            }
        }
    } catch (const std::exception& ex) {
        std::cerr << "Input Error: " << ex.what() << "\n";
        return 1;
    }

    const std::size_t n_nodes = gammas.size();
    const std::vector<long long> primes = spectral::generate_primes(n_nodes);

    std::vector<spectral::SpectralNode> riemann_nodes;
    std::vector<spectral::SpectralNode> baseline_nodes;
    riemann_nodes.reserve(n_nodes);
    baseline_nodes.reserve(n_nodes);

    std::mt19937_64 base_rng(cfg.seed);
    std::uniform_real_distribution<double> phase_dist(0.0, 2.0 * spectral::PI);
    spectral::PhaseAccumulatorState split_state;

    for (std::size_t i = 0; i < n_nodes; ++i) {
        long long p = primes[i];
        double gamma = gammas[i];

        // 1. Riemann Zeta Phase Projection + Q Lift
        double r_theta = std::fmod(gamma * std::log(static_cast<double>(p)), 2.0 * spectral::PI);
        if (r_theta < 0.0) r_theta += 2.0 * spectral::PI;
        long long r_d_init = std::llround(static_cast<double>(p) * 0.5 * std::cos(r_theta));
        auto r_node = spectral::resolve_from_initial_d(i, gamma, p, r_theta, r_d_init, cfg.max_steps);

        // 2. Intrinsic Prime-Split Brahmagupta-Fibonacci Phase Accumulator (BigRational S^1(Q))
        r_node.split_phase = spectral::advance_split_phase(split_state, p, r_theta);

        // 3. Seeded Random Phase Baseline Projection + Q Lift
        double b_theta = phase_dist(base_rng);
        long long b_d_init = std::llround(static_cast<double>(p) * 0.5 * std::cos(b_theta));
        auto b_node = spectral::resolve_from_initial_d(i, gamma, p, b_theta, b_d_init, cfg.max_steps);

        if (!verify_node_invariants(r_node, "Riemann") ||
            !verify_node_invariants(b_node, "Baseline") ||
            !spectral::verify_split_phase_invariants(p, r_node.split_phase)) {
            return 2;
        }

        riemann_nodes.push_back(std::move(r_node));
        baseline_nodes.push_back(std::move(b_node));
    }

    if (!spectral::export_side_by_side_csv(cfg.output_csv, riemann_nodes, baseline_nodes)) {
        std::cerr << "Error: Failed to write CSV to " << cfg.output_csv << "\n";
        return 1;
    }

    if (!cfg.quiet) {
        std::cout << "Idx | Prime | Gamma    | R_D0 -> R_Dc (St) | R_Q_Det | Split_Theta | Cum_Wrap | Phase_Diff | S1(Q) Digits\n"
                  << "----+-------+----------+-------------------+---------+-------------+----------+------------+-------------\n";
        for (std::size_t i = 0; i < n_nodes; ++i) {
            const auto& r = riemann_nodes[i];
            const auto& sp = r.split_phase;
            std::cout << std::setw(3) << r.index << " | "
                      << std::setw(5) << r.prime << " | "
                      << std::fixed << std::setprecision(3) << std::setw(8) << r.riemann_zero_gamma << " | "
                      << std::setw(4) << r.d_initial << " -> " << std::setw(4) << r.d_converged
                      << " (" << std::setw(2) << r.search_steps << ") | "
                      << std::setw(7) << r.rat.norm_det.str() << " | "
                      << std::setprecision(4) << std::setw(11) << sp.step_theta << " | "
                      << std::setw(8) << sp.cumulative_wrapped << " | "
                      << std::showpos << std::setw(10) << sp.riemann_phase_diff << std::noshowpos << " | "
                      << std::setw(12) << sp.q_re.compact_str() << "\n";
        }
    }

    spectral::MonteCarloSummary mc = spectral::run_monte_carlo_comparison(
        riemann_nodes, cfg.max_steps, cfg.mc_trials, cfg.seed + 1ULL);

    const auto& final_sp = riemann_nodes.back().split_phase;
    std::cout << "\n=== Coupled Spectral Lattice & Brahmagupta-Fibonacci Phase Report ===\n"
              << "Nodes Evaluated:          " << mc.node_count << " (p_max = " << primes.back() << ")\n"
              << "Max Search Radius:        " << cfg.max_steps << "\n"
              << "Riemann Exact Matches:    " << mc.riemann_exact_matches << " / " << mc.node_count
              << " (" << std::fixed << std::setprecision(2)
              << (100.0 * mc.riemann_exact_matches / mc.node_count) << "%)\n"
              << "Riemann Mean Steps:       " << mc.riemann_mean_steps << "\n"
              << "Split Phase Unwrapped:    " << std::setprecision(6) << final_sp.cumulative_total << " rad\n"
              << "Split Phase Wrapped:      " << final_sp.cumulative_wrapped << " rad ("
              << (final_sp.cumulative_wrapped * 180.0 / spectral::PI) << " deg)\n"
              << "Final S^1(Q) Precision:   " << final_sp.q_re.den.str().size() << " decimal digits (Exact Unit Norm: VERIFIED)\n";

    if (cfg.mc_trials > 0) {
        std::cout << "Monte Carlo Trials:       " << mc.trials << "\n"
                  << "Baseline Mean Exact:      " << mc.baseline_mean_exact << " +/- " << mc.baseline_std_exact
                  << " (" << (100.0 * mc.baseline_mean_exact / mc.node_count) << "%)\n"
                  << "Baseline Mean Steps:      " << mc.baseline_mean_steps << "\n"
                  << "Exact Match Z-Score:      " << std::showpos << mc.z_score_exact << std::noshowpos << "\n"
                  << "Empirical One-Sided p:    " << std::setprecision(4) << mc.empirical_p_value << "\n";
    }
    std::cout << "Side-by-Side CSV Saved:   " << cfg.output_csv << "\n";

    if (cfg.verify_gue) {
        spectral::verify_gue_semicircle_unfolding(cfg.gue_size, cfg.seed);
    }

    return 0;
}
