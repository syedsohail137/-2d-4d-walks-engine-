#include "manifold.h"

bool verify_norm_multiplicativity(int16_t a1, int16_t b1, int16_t c1, int16_t d1,
                                  int16_t a2, int16_t b2, int16_t c2, int16_t d2, 
                                  uint64_t p1, uint64_t p2) {
    int64_t target_norm = static_cast<int64_t>(p1) * static_cast<int64_t>(p2);
    int64_t rs_a = (a1*a2 - b1*b2 - c1*c2 - d1*d2);
    int64_t rs_b = (a1*b2 + b1*a2 + c1*d2 - d1*c2);
    int64_t rs_c = (a1*c2 - b1*d2 + c1*a2 + d1*b2);
    int64_t rs_d = (a1*d2 + b1*c2 - c1*b2 + d1*a2);
    int64_t computed_norm = (rs_a*rs_a) + (rs_b*rs_b) + (rs_c*rs_c) + (rs_d*rs_d);
    return computed_norm == target_norm;
}
