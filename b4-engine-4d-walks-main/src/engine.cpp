#include "engine.hpp"
#include "math_bounds.hpp"
#include <omp.h>
#include <random>
#include <cmath>

SimulationResult execute_simulation(size_t total_walks, size_t max_steps) {
    double start_time = omp_get_wtime();
    
    double total_squared_distance_accumulator = 0.0;
    size_t total_absorbed_walks_accumulator = 0;

    // Parallel execution block with local variable arithmetic reductions
    #pragma omp parallel reduction(+:total_squared_distance_accumulator, total_absorbed_walks_accumulator)
    {
        int thread_id = omp_get_thread_num();
        
        // Thread-isolated PCG/Mersenne-Twister engine avoids atomic locks and heap false sharing
        std::mt19937_64 prng(1337 + thread_id); 
        
        // Distribution bounds mapping to 4 discrete dimensions and binary directions
        std::uniform_int_distribution<size_t> coord_dist(0, 3);
        std::uniform_int_distribution<int> direction_dist(0, 1);

        // Uniform workload allocation across threads with non-overlapping iteration ranges
        #pragma omp for schedule(static)
        for (size_t i = 0; i < total_walks; ++i) {
            Position4D walk_pos;
            bool is_absorbed = false;

            for (size_t step = 0; step < max_steps; ++step) {
                size_t TargetDimension = coord_dist(prng);
                
                // Map binary uniform output {0, 1} strictly to spatial steps {-1, +1}
                int StepDirection = direction_dist(prng) * 2 - 1; 
                
                walk_pos.coords[TargetDimension] += StepDirection;

                // Continuous Robbins-Siegmund energy containment validation check
                if (evaluate_absorption_criteria(walk_pos)) {
                    is_absorbed = true;
                    break; 
                }
            }

            // Extract terminal Euclidean metrics
            double individual_sq_dist = 0.0;
            for (int d = 0; d < 4; ++d) {
                individual_sq_dist += static_cast<double>(walk_pos.coords[d]) * walk_pos.coords[d];
            }

            total_squared_distance_accumulator += individual_sq_dist;
            if (is_absorbed) {
                total_absorbed_walks_accumulator++;
            }
        }
    }

    double end_time = omp_get_wtime();
    
    // Compute core structural tracking laws: <||X(t)||_2> = alpha * sqrt(t)
    double mean_squared_distance = total_squared_distance_accumulator / static_cast<double>(total_walks);
    double rms_displacement = std::sqrt(mean_squared_distance);
    double empirical_alpha = rms_displacement / std::sqrt(static_cast<double>(max_steps));

    SimulationResult metrics;
    metrics.root_mean_square_displacement = rms_displacement;
    metrics.empirical_drift_coefficient = empirical_alpha;
    metrics.absorption_count = total_absorbed_walks_accumulator;
    metrics.total_execution_time = end_time - start_time;

    return metrics;
}
