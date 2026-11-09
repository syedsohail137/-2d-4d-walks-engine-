#ifndef B4_ENGINE_HPP
#define B4_ENGINE_HPP

#include <cstdint>
#include <vector>

// Explicitly aligned to 32 bytes to ensure it cannot span across L1/L2 cache line boundaries
struct alignas(32) Position4D {
    int32_t coords[4] = {0, 0, 0, 0};
};

// Global simulation results matrix payload
struct SimulationResult {
    double root_mean_square_displacement;
    double empirical_drift_coefficient;
    size_t absorption_count;
    double total_execution_time;
};

// Core high-performance simulation execution interface
SimulationResult execute_simulation(size_t total_walks, size_t max_steps);

#endif // B4_ENGINE_HPP
