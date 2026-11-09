#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <cstdint>
#include <vector>

struct alignas(8) PrimeNode4D {
    uint64_t prime;
    int16_t a, b, c, d;
    float phase_r;
    float phase_i;
};

bool decompose_four_square_optimized(uint64_t p, int16_t& a, int16_t& b, int16_t& c, int16_t& d);
bool verify_norm_multiplicativity(int16_t a1, int16_t b1, int16_t c1, int16_t d1,
                                  int16_t a2, int16_t b2, int16_t c2, int16_t d2, 
                                  uint64_t p1, uint64_t p2);

#endif
