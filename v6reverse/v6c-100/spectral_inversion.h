#ifndef SPECTRAL_INVERSION_H
#define SPECTRAL_INVERSION_H

#include <cstdint>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

namespace spectral {

constexpr double PI = 3.14159265358979323846;

struct Rational {
    long long num = 0;
    long long den = 1;

    Rational(long long n = 0, long long d = 1) {
        if (d == 0) { num = 0; den = 1; return; }
        if (d < 0) { n = -n; d = -d; }
        long long g = std::gcd(n < 0 ? -n : n, d);
        num = (g > 0) ? (n / g) : 0;
        den = (g > 0) ? (d / g) : 1;
    }

    double to_double() const {
        return static_cast<double>(num) / static_cast<double>(den);
    }

    std::string str() const {
        if (den == 1) return std::to_string(num);
        return std::to_string(num) + "/" + std::to_string(den);
    }
};

inline std::ostream& operator<<(std::ostream& os, const Rational& r) {
    return os << r.str();
}

struct QuaternionLattice {
    long long a = 0;
    long long b = 0;
    long long c = 0;
    long long d = 0;
};

struct RationalProjection {
    Rational norm_det;   // D_converged / p in Q
    Rational norm_sq;    // (a^2 + b^2 + c^2 + d^2) / p^2 = 1/p in Q
    Rational stereo_x;   // a / (p - d) in Q
    Rational stereo_y;   // b / (p - d) in Q
    Rational stereo_z;   // c / (p - d) in Q
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
    RationalProjection rat;
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
    std::size_t count = 100; // Defaults to all 100 embedded zeros (use --count 50 for base case)
    long long max_steps = 40;
    int mc_trials = 500;
    uint64_t seed = 42ULL;
    std::string gammas_file;
    std::string output_csv = "spectral_comparison_100.csv";
    bool quiet = false;
    bool verify_gue = false;
    int gue_size = 200;
};

// Verified Table: First 100 non-trivial zeros of the Riemann Zeta function (Im(s) > 0)
inline const std::vector<double> DEFAULT_RIEMANN_GAMMAS = {
    // Zeros 1 - 50 (Base Case)
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
    141.123707404021124, 143.111845807620632,
    // Zeros 51 - 100 (Extended Scaling to Q)
    146.000982486765519, 147.422765342559602, 150.053520420784880,
    150.925257612241467, 153.024693811198896, 156.112909294237868,
    157.597591817594060, 158.849988171420499, 161.188964137596028,
    163.030709687181987, 165.537069187900419, 167.184439978174513,
    169.094515415568821, 169.911976479411699, 173.411536519591553,
    174.754191523365726, 176.441434297710419, 178.377407776099977,
    179.916484020256996, 182.207078484366462, 184.874467848387509,
    185.598783677707471, 187.228922583501852, 189.416158656016937,
    192.026656360713787, 193.079726603845704, 195.265396679529235,
    196.876481840958317, 198.015309676251912, 201.264751943703789,
    202.493594514140534, 204.189671803104554, 205.394697202163286,
    207.906258887806210, 209.576509716856260, 211.690862595365308,
    213.347919359712666, 214.547044783491423, 216.169538508263700,
    219.067596349021379, 220.714918839314003, 221.430705554693339,
    224.007000254604335, 224.983324669582288, 227.421444279679291,
    229.337413305525348, 231.250188700499165, 231.987235253180249,
    233.693404178908301, 236.524229665816206
};

} // namespace spectral

#endif // SPECTRAL_INVERSION_H
