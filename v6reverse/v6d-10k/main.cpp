#include "io_and_montecarlo.hpp"
#include "lattice_solver.hpp"
#include "spectral_inversion.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "Options:\n"
              << "  --gammas-file <path>      Load Riemann zeros from external LMFDB/Odlyzko file\n"
              << "  --projection <mode>       Phase projection geometry: cosine | linear (default: cosine)\n"
              << "  --count <N>               Number of prime-zero pairs to evaluate (default: all loaded)\n"
              << "  --max-steps <N>           Maximum determinant search radius (default: 100)\n"
              << "  --mc-trials <N>           Monte Carlo random phase trials for p-value (default: 1000)\n"
              << "  --threads <N>             Number of OpenMP threads (default: 0 = auto)\n"
              << "  --seed <N>                64-bit RNG seed for baseline & Monte Carlo (default: 42)\n"
              << "  --output <path>           Side-by-side comparison CSV path (default: spectral_comparison_large.csv)\n"
              << "  --quiet                   Suppress per-row stdout preview\n"
              << "  --help                    Display this help message\n";
}

bool parse_cli(int argc, char* argv[], spectral::RunConfig& cfg) {
    const char* prog_name = (argc > 0 && argv != nullptr && argv[0] != nullptr) ? argv[0] : "spectral_inversion";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(prog_name);
            std::exit(0);
        } else if (arg == "--quiet") {
            cfg.quiet = true;
        } else if (arg == "--gammas-file" && i + 1 < argc) {
            cfg.gammas_file = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            cfg.output_csv = argv[++i];
        } else if (arg == "--projection" && i + 1 < argc) {
            std::string mode = argv[++i];
            if (mode == "cosine") {
                cfg.projection = spectral::ProjectionMode::Cosine;
            } else if (mode == "linear") {
                cfg.projection = spectral::ProjectionMode::Linear;
            } else {
                std::cerr << "Error: --projection must be 'cosine' or 'linear'.\n";
                return false;
            }
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
        } else if (arg == "--threads" && i + 1 < argc) {
            cfg.threads = std::stoi(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            cfg.seed = std::stoull(argv[++i]);
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            print_usage(prog_name);
            return false;
        }
    }
    return true;
}

bool verify_node_invariants(const spectral::SpectralNode& n, const char* label) {
    if (n.state == "No Representation") return true;
    long long norm_sq = n.q.a * n.q.a + n.q.b * n.q.b + n.q.c * n.q.c + n.q.d * n.q.d;
    long long det = n.q.a * n.q.d - n.q.b * n.q.c;
    if (norm_sq != n.prime || det != n.d_converged) {
        std::cerr << "Algebraic invariant failure [" << label << "] at Index " << n.index << "\n";
        return false;
    }
    return true;
}

void print_scale_tiers(const std::vector<spectral::SpectralNode>& riemann_nodes,
                       const std::vector<spectral::SpectralNode>& baseline_nodes) {
    const std::vector<std::pair<std::size_t, std::size_t>> tiers = {
        {0, 100}, {100, 1000}, {1000, 5000}, {5000, 10000}, {10000, riemann_nodes.size()}
    };

    std::cout << "\n--- Magnitude Scale-Tier Breakdown ---\n"
              << "Index Range     | Prime Range        | Riemann Exact | Single-Seed Base Exact | Mean Steps (R / B)\n"
              << "----------------+--------------------+---------------+------------------------+-------------------\n";

    for (const auto& tier : tiers) {
        std::size_t start = tier.first;
        std::size_t end = std::min(tier.second, riemann_nodes.size());
        if (start >= end) continue;

        std::size_t count = end - start;
        std::size_t r_exact = 0, b_exact = 0;
        long long r_steps = 0, b_steps = 0;

        for (std::size_t i = start; i < end; ++i) {
            if (riemann_nodes[i].state == "Exact Match") ++r_exact;
            if (baseline_nodes[i].state == "Exact Match") ++b_exact;
            r_steps += riemann_nodes[i].search_steps;
            b_steps += baseline_nodes[i].search_steps;
        }

        std::cout << std::setw(6) << (start + 1) << ".." << std::setw(7) << std::left << end << std::right << " | "
                  << std::setw(7) << riemann_nodes[start].prime << ".." << std::setw(9) << std::left << riemann_nodes[end - 1].prime << std::right << " | "
                  << std::setw(6) << std::fixed << std::setprecision(2) << (100.0 * r_exact / count) << "%       | "
                  << std::setw(6) << (100.0 * b_exact / count) << "%                | "
                  << std::setprecision(2) << (static_cast<double>(r_steps) / count)
                  << " / " << (static_cast<double>(b_steps) / count) << "\n";
    }
}

} // namespace

int main(int argc, char* argv[]) {
    spectral::RunConfig cfg;
    if (!parse_cli(argc, argv, cfg)) return 1;

    int active_threads = 1;
#ifdef _OPENMP
    if (cfg.threads > 0) {
        omp_set_num_threads(cfg.threads);
    }
    active_threads = omp_get_max_threads();
#endif

    auto t_start = std::chrono::high_resolution_clock::now();

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

    std::vector<spectral::SpectralNode> riemann_nodes(n_nodes);
    std::vector<spectral::SpectralNode> baseline_nodes(n_nodes);

    std::vector<double> base_phases(n_nodes);
    {
        std::mt19937_64 base_rng(cfg.seed);
        std::uniform_real_distribution<double> phase_dist(0.0, 2.0 * spectral::PI);
        for (std::size_t i = 0; i < n_nodes; ++i) {
            base_phases[i] = phase_dist(base_rng);
        }
    }

    #pragma omp parallel
    {
        std::vector<std::pair<long long, long long>> buf1;
        std::vector<std::pair<long long, long long>> buf2;
        buf1.reserve(32);
        buf2.reserve(32);

        #pragma omp for schedule(dynamic, 32)
        for (std::size_t i = 0; i < n_nodes; ++i) {
            long long p = primes[i];
            double gamma = gammas[i];

            double r_theta = std::fmod(gamma * std::log(static_cast<double>(p)), 2.0 * spectral::PI);
            if (r_theta < 0.0) r_theta += 2.0 * spectral::PI;
            long long r_d_init = spectral::project_phase_to_determinant(p, r_theta, cfg.projection);
            riemann_nodes[i] = spectral::resolve_from_initial_d(
                i, gamma, p, r_theta, r_d_init, cfg.max_steps, buf1, buf2);

            double b_theta = base_phases[i];
            long long b_d_init = spectral::project_phase_to_determinant(p, b_theta, cfg.projection);
            baseline_nodes[i] = spectral::resolve_from_initial_d(
                i, gamma, p, b_theta, b_d_init, cfg.max_steps, buf1, buf2);
        }
    }

    for (std::size_t i = 0; i < n_nodes; ++i) {
        if (!verify_node_invariants(riemann_nodes[i], "Riemann") ||
            !verify_node_invariants(baseline_nodes[i], "Baseline")) {
            return 2;
        }
    }

    if (!spectral::export_side_by_side_csv(cfg.output_csv, riemann_nodes, baseline_nodes)) {
        std::cerr << "Error: Failed to write CSV to " << cfg.output_csv << "\n";
        return 1;
    }

    if (!cfg.quiet) {
        std::size_t preview = std::min<std::size_t>(n_nodes, 30);
        std::cout << "Previewing first " << preview << " of " << n_nodes << " rows:\n"
                  << "Idx | Prime | Gamma    | R_D0 -> R_Dc (Steps) | B_D0 -> B_Dc (Steps)\n"
                  << "----+-------+----------+----------------------+---------------------\n";
        for (std::size_t i = 0; i < preview; ++i) {
            const auto& r = riemann_nodes[i];
            const auto& b = baseline_nodes[i];
            std::cout << std::setw(3) << r.index << " | "
                      << std::setw(5) << r.prime << " | "
                      << std::fixed << std::setprecision(3) << std::setw(8) << r.riemann_zero_gamma << " | "
                      << std::setw(4) << r.d_initial << " -> " << std::setw(4) << r.d_converged
                      << " (" << std::setw(2) << r.search_steps << ")   | "
                      << std::setw(4) << b.d_initial << " -> " << std::setw(4) << b.d_converged
                      << " (" << std::setw(2) << b.search_steps << ")\n";
        }
    }

    spectral::MonteCarloSummary mc = spectral::run_monte_carlo_comparison(
        riemann_nodes, cfg.max_steps, cfg.mc_trials, cfg.seed + 1ULL, cfg.projection);

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_sec = std::chrono::duration<double>(t_end - t_start).count();

    print_scale_tiers(riemann_nodes, baseline_nodes);

    std::cout << "\n=== OpenMP Spectral Lattice Projection & Monte Carlo Report ===\n"
              << "Projection Geometry:      "
              << (cfg.projection == spectral::ProjectionMode::Linear ? "Linear (Area-Preserving)" : "Cosine") << "\n"
              << "OpenMP Active Threads:    " << active_threads << "\n"
              << "Wall-Clock Elapsed:       " << std::fixed << std::setprecision(3) << elapsed_sec << " s\n"
              << "Nodes Evaluated:          " << mc.node_count << " (p_max = " << primes.back() << ")\n"
              << "Max Search Radius:        " << cfg.max_steps << "\n"
              << "Riemann Exact Matches:    " << mc.riemann_exact_matches << " / " << mc.node_count
              << " (" << std::setprecision(2)
              << (100.0 * mc.riemann_exact_matches / mc.node_count) << "%)\n"
              << "Riemann Mean Steps:       " << mc.riemann_mean_steps << "\n";

    if (cfg.mc_trials > 0) {
        std::cout << "Monte Carlo Trials:       " << mc.trials << "\n"
                  << "Baseline Mean Exact:      " << mc.baseline_mean_exact << " +/- " << mc.baseline_std_exact
                  << " (" << (100.0 * mc.baseline_mean_exact / mc.node_count) << "%)\n"
                  << "Baseline Mean Steps:      " << mc.baseline_mean_steps << "\n"
                  << "Exact Match Z-Score:      " << std::showpos << mc.z_score_exact << std::noshowpos << "\n"
                  << "Empirical One-Sided p:    " << std::setprecision(4) << mc.empirical_p_value << "\n";
    }
    std::cout << "Side-by-Side CSV Saved:   " << cfg.output_csv << "\n";

    return 0;
}
