#ifndef RATIONAL_PROJECTION_HPP
#define RATIONAL_PROJECTION_HPP

#include "spectral_inversion.h"

namespace spectral {

// Lifts a solved Lipschitz quaternion q = (a, b, c, d) with a^2+b^2+c^2+d^2 = p
// onto exact Q coordinates:
//   1. Normalized rational determinant: D_conv / p in Q
//   2. Unit rational quaternion U = q^2 / p in S^3(Q) with exact norm ||U||^2 = 1
//   3. Stereographic coordinates (a/(p-d), b/(p-d), c/(p-d)) in Q^3
inline RationalProjection lift_to_rationals(long long p, long long d_conv, const QuaternionLattice& q) {
    RationalProjection rp;
    if (p <= 0) return rp;

    // 1. Normalized rational determinant D/p in [-1/2, 1/2]
    rp.norm_det = Rational(d_conv, p);

    // 2. Exact unit norm in S^3(Q) via rational quaternion squaring U = q^2 / p:
    //    U = ((a^2 - b^2 - c^2 - d^2)/p, 2ab/p, 2ac/p, 2ad/p)
    long long u0 = q.a * q.a - q.b * q.b - q.c * q.c - q.d * q.d;
    long long u1 = 2 * q.a * q.b;
    long long u2 = 2 * q.a * q.c;
    long long u3 = 2 * q.a * q.d;
    long long unit_norm_num = u0 * u0 + u1 * u1 + u2 * u2 + u3 * u3;
    rp.norm_sq = Rational(unit_norm_num, p * p); // Always reduces to exact 1/1 in Q

    // 3. Rational stereographic coordinates in Q^3
    long long denom = p - q.d;
    rp.stereo_x = Rational(q.a, denom);
    rp.stereo_y = Rational(q.b, denom);
    rp.stereo_z = Rational(q.c, denom);

    return rp;
}

// Verifies exact rational field invariants in Q without floating-point rounding
inline bool verify_rational_invariants(const SpectralNode& n) {
    if (n.state == "No Representation") return true;

    // Check 1: Unit rational quaternion norm ||q^2 / p||^2 == 1/1 in Q
    if (n.rat.norm_sq.num != 1 || n.rat.norm_sq.den != 1) {
        return false;
    }

    // Check 2: Exact stereographic identity on S^3(sqrt(p)):
    //          (a/(p-d))^2 + (b/(p-d))^2 + (c/(p-d))^2 == (p - d^2) / (p - d)^2
    long long lhs_num = n.q.a * n.q.a + n.q.b * n.q.b + n.q.c * n.q.c;
    long long rhs_num = n.prime - (n.q.d * n.q.d);

    return lhs_num == rhs_num;
}

} // namespace spectral

#endif // RATIONAL_PROJECTION_HPP
