#include "manifold.h"
#include <cmath>
#include <algorithm>

bool decompose_four_square_optimized(uint64_t p, int16_t& a, int16_t& b, int16_t& c, int16_t& d) {
    int64_t max_a = static_cast<int64_t>(std::sqrt(p));
    for (int64_t ta = max_a; ta >= 0; --ta) {
        int64_t rem_a = p - (ta * ta); if (rem_a < 0) continue;
        int64_t max_b = static_cast<int64_t>(std::sqrt(rem_a));
        for (int64_t tb = max_b; tb >= 0; --tb) {
            int64_t rem_b = rem_a - (tb * tb); if (rem_b < 0) continue;
            int64_t max_c = static_cast<int64_t>(std::sqrt(rem_b));
            for (int64_t tc = max_c; tc >= 0; --tc) {
                int64_t rem_c = rem_b - (tc * tc); if (rem_c < 0) continue;
                int64_t td = static_cast<int64_t>(std::round(std::sqrt(rem_c)));
                if (ta * ta + tb * tb + tc * tc + td * td == p) {
                    int16_t coords[4] = { static_cast<int16_t>(ta), static_cast<int16_t>(tb), 
                                          static_cast<int16_t>(tc), static_cast<int16_t>(td) };
                    std::sort(coords, coords + 4, std::greater<int16_t>());
                    if (p % 2 != 0 && coords[0] % 2 == 0) {
                        for (int i = 1; i < 4; ++i) {
                            if (coords[i] % 2 != 0) {
                                std::swap(coords[0], coords[i]);
                                break;
                            }
                        }
                    }
                    a = coords[0]; b = coords[1]; c = coords[2]; d = coords[3];
                    return true;
                }
            }
        }
    }
    return false;
}
