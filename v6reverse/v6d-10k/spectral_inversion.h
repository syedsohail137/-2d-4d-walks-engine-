#ifndef SPECTRAL_INVERSION_H
#define SPECTRAL_INVERSION_H

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace spectral {

constexpr double PI = 3.14159265358979323846;

enum class ProjectionMode {
    Cosine,
    Linear
};

struct QuaternionLattice {
    long long a = 0;
    long long b = 0;
    long long c = 0;
    long long d = 0;
};

struct SpectralNode {
    int index = 0;
    double riemann_zero_gamma = 0.0;
    long long prime = 0;
    double phase_angle_rad = 0.0;
    long long d_initial = 0;
    long long d_converged = 0;
    int search_steps = 0;
    long long adjustment = 0;
    std::string state;
    QuaternionLattice q;
};

struct MonteCarloSummary {
    int trials = 0;
    std::size_t node_count = 0;
    std::size_t riemann_exact_matches = 0;
    double riemann_mean_steps = 0.0;
    double baseline_mean_exact = 0.0;
    double baseline_std_exact = 0.0;
    double baseline_mean_steps = 0.0;
    double z_score_exact = 0.0;
    double empirical_p_value = 0.0;
};

struct RunConfig {
    std::size_t count = 0;
    long long max_steps = 100;
    int mc_trials = 1000;
    int threads = 0;
    uint64_t seed = 42ULL;
    ProjectionMode projection = ProjectionMode::Cosine;
    std::string gammas_file;
    std::string output_csv = "spectral_comparison_large.csv";
    bool quiet = false;
};

// Maps phase theta in [0, 2*PI) to an initial determinant in [-p/2, p/2]
inline long long project_phase_to_determinant(long long p, double theta, ProjectionMode mode) {
    if (mode == ProjectionMode::Linear) {
        double normalized = (theta / (2.0 * PI)) - 0.5;
        return std::llround(static_cast<double>(p) * normalized);
    }
    return std::llround(static_cast<double>(p) * 0.5 * std::cos(theta));
}

inline const std::vector<double> DEFAULT_RIEMANN_GAMMAS = {
    14.134725141734693, 21.022039638771554, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302, 52.970321477714460, 56.446247697063394,
    59.347044002602353, 60.831778524609878, 65.112544048081606,
    67.079810529494173, 69.546401711173970, 72.067157674481270,
    75.704690699083430, 77.144840068874000, 79.337375020249000,
    82.910380854086000, 84.735492980500050, 87.425274613957630,
    88.809111207634460, 92.491899270558000, 94.651344040597000,
    95.870634228245000, 98.831194287092000, 101.317851005796000,
    103.725538040478334, 105.446623052326091, 107.168611184276408,
    111.029535543169674, 111.874659176992637, 114.320220915452713,
    116.226680320857554, 118.790782865976217, 121.370125002420646,
    122.946829293552588, 124.256818554345767, 127.516683879596495,
    129.578704199956051, 131.087688530932657, 133.497737202997586,
    134.756509753373871, 138.116042054533443, 139.736208952121389,
    141.123707404021124, 143.111845807620632
};

} // namespace spectral

#endif // SPECTRAL_INVERSION_H
