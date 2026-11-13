#ifndef GUE_UNFOLDING_HPP
#define GUE_UNFOLDING_HPP

#include "spectral_inversion.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

namespace spectral {

// Householder tridiagonalization + Implicit QL eigensolver for complex Hermitian matrices
inline std::vector<double> compute_hermitian_eigenvalues(
    std::vector<std::vector<std::complex<double>>> H) {

    int n = static_cast<int>(H.size());
    std::vector<double> d(n, 0.0);
    std::vector<double> e(n, 0.0);

    for (int i = n - 1; i >= 1; --i) {
        double scale = 0.0, h = 0.0;
        for (int k = 0; k < i; ++k) {
            scale += std::abs(H[i][k].real()) + std::abs(H[i][k].imag());
        }
        if (scale == 0.0) {
            e[i] = 0.0;
        } else {
            for (int k = 0; k < i; ++k) {
                H[i][k] /= scale;
                h += std::norm(H[i][k]);
            }
            std::complex<double> f = H[i][i - 1];
            double abs_f = std::abs(f);
            std::complex<double> g = (abs_f > 0.0)
                ? (-std::sqrt(h) * (f / abs_f))
                : std::complex<double>(-std::sqrt(h), 0.0);
            e[i] = scale * std::abs(g);
            h -= abs_f * std::sqrt(h);
            H[i][i - 1] = f - g;

            std::vector<std::complex<double>> u(i), p(i, {0.0, 0.0});
            for (int j = 0; j < i; ++j) u[j] = H[i][j] / h;

            std::complex<double> u_dot_p(0.0, 0.0);
            for (int j = 0; j < i; ++j) {
                for (int k = 0; k <= j; ++k) p[j] += H[j][k] * u[k];
                for (int k = j + 1; k < i; ++k) p[j] += std::conj(H[k][j]) * u[k];
                u_dot_p += std::conj(u[j]) * p[j];
            }
            std::complex<double> K = 0.5 * u_dot_p;
            for (int j = 0; j < i; ++j) p[j] -= K * u[j];

            for (int j = 0; j < i; ++j) {
                for (int k = 0; k <= j; ++k) {
                    H[j][k] -= (u[j] * std::conj(p[k]) + p[j] * std::conj(u[k]));
                }
            }
        }
        d[i] = H[i][i].real();
    }
    d[0] = H[0][0].real();

    for (int i = 0; i < n - 1; ++i) e[i] = e[i + 1];
    e[n - 1] = 0.0;

    for (int l = 0; l < n; ++l) {
        int iter = 0, m;
        do {
            for (m = l; m < n - 1; ++m) {
                double dd = std::abs(d[m]) + std::abs(d[m + 1]);
                if (std::abs(e[m]) + dd == dd) break;
            }
            if (m != l) {
                if (++iter > 60) break;
                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = std::hypot(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + std::copysign(r, g));
                double s = 1.0, c = 1.0, p = 0.0;
                int i;
                for (i = m - 1; i >= l; --i) {
                    double f = s * e[i], b = c * e[i];
                    r = std::hypot(f, g);
                    e[i + 1] = r;
                    if (r == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
                    s = f / r; c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p; e[l] = g; e[m] = 0.0;
            }
        } while (m != l);
    }

    std::sort(d.begin(), d.end());
    return d;
}

// Corrected Module 3: Unfolds GUE spectrum via exact Wigner Semicircle CDF F(lambda)
// Restores unit local density and yields ~40% peak concentration in 0.75 < s < 1.25
inline void verify_gue_semicircle_unfolding(int matrix_size, uint64_t seed) {
    std::cout << "\n=== Module 3: GUE Semicircle CDF Unfolding Verification (N = " << matrix_size << ") ===\n";

    std::mt19937_64 prng(seed);
    std::normal_distribution<double> gaussian(0.0, 1.0);

    std::vector<std::vector<std::complex<double>>> H(
        matrix_size, std::vector<std::complex<double>>(matrix_size));

    for (int i = 0; i < matrix_size; ++i) {
        for (int j = i; j < matrix_size; ++j) {
            if (i == j) {
                H[i][j] = std::complex<double>(gaussian(prng), 0.0);
            } else {
                double re = gaussian(prng) / std::sqrt(2.0);
                double im = gaussian(prng) / std::sqrt(2.0);
                H[i][j] = std::complex<double>(re, im);
                H[j][i] = std::conj(H[i][j]);
            }
        }
    }

    std::vector<double> eigenvalues = compute_hermitian_eigenvalues(H);

    // Wigner semicircle spectral radius R = 2 * sqrt(N)
    const double R = 2.0 * std::sqrt(static_cast<double>(matrix_size));
    std::vector<double> unfolded(matrix_size);

    for (int i = 0; i < matrix_size; ++i) {
        double x = std::clamp(eigenvalues[i], -R, R);
        double cdf = 0.5
                   + (x * std::sqrt(std::max(0.0, R * R - x * x))) / (2.0 * PI * R * R)
                   + std::asin(x / R) / PI;
        unfolded[i] = static_cast<double>(matrix_size) * cdf;
    }

    const int trim = static_cast<int>(matrix_size * 0.15);
    std::vector<double> spacings;
    spacings.reserve(matrix_size - 2 * trim - 1);

    double sum_s = 0.0;
    for (int i = trim; i < matrix_size - trim - 1; ++i) {
        double s = unfolded[i + 1] - unfolded[i];
        spacings.push_back(s);
        sum_s += s;
    }

    const double mean_s = sum_s / static_cast<double>(spacings.size());
    int close_count = 0;
    int peak_count = 0;

    for (double& s : spacings) {
        s /= mean_s;
        if (s < 0.15) ++close_count;
        if (s > 0.75 && s < 1.25) ++peak_count;
    }

    const double total = static_cast<double>(spacings.size());
    std::cout << "Unfolded Bulk Spacings:   " << spacings.size() << "\n"
              << "Close Levels (s < 0.15):  " << close_count
              << " (" << std::fixed << std::setprecision(1) << (100.0 * close_count / total) << "% | Theory: ~0.4%)\n"
              << "Peak Region (0.75..1.25): " << peak_count
              << " (" << std::fixed << std::setprecision(1) << (100.0 * peak_count / total) << "% | Theory: ~40.2%)\n";
}

} // namespace spectral

#endif // GUE_UNFOLDING_HPP
