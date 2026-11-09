#include "engine.hpp"
#include <iostream>
#include <iomanip>

int main() {
    // Execution parameters locked to theoretical validation requirements
    const size_t TOTAL_WALKS = 1000000; 
    const size_t MAX_STEPS   = 5000;   

    std::cout << "========================================================\n";
    std::cout << " INITIALIZING B4_ENGINE: 4D MULTI-THREADED WALKS ENGINE \n";
    std::cout << "========================================================\n";
    std::cout << "Executing: " << TOTAL_WALKS << " parallel independent paths.\n";
    std::cout << "Max Step Boundary Depth per path: " << MAX_STEPS << " steps.\n\n";

    // Run high performance calculation loop
    SimulationResult output = execute_simulation(TOTAL_WALKS, MAX_STEPS);

    // Print empirical metrics compared against theoretical benchmarks
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "---------------- SIMULATION METRICS --------------------\n";
    std::cout << "Execution Multi-threading Duration : " << output.total_execution_time << " seconds\n";
    std::cout << "Total Boundary Absorbed Walks      : " << output.absorption_count << " / " << TOTAL_WALKS << "\n";
    std::cout << "Empirical Root-Mean-Square <||X||> : " << output.root_mean_square_displacement << "\n";
    std::cout << "Calculated Scaling Exponent Alpha  : " << output.empirical_drift_coefficient << "\n";
    std::cout << "---------------- VALIDATION CHECKS ---------------------\n";
    
    // Cross check empirical scaling coefficient against 0.8 * sqrt(t) target bound
    std::cout << "Theoretical Target Exponent Bound  : ~0.80000\n";
    double variance = std::abs(output.empirical_drift_coefficient - 0.8);
    std::cout << "Absolute Coefficient Delta Variance: " << variance << "\n";
    
    if (variance <= 0.05) {
        std::cout << "STATUS: SUCCESS [Empirical metrics fit asymptotic limits]\n";
    } else {
        std::cout << "STATUS: UNVERIFIED [Check continuous absorption boundaries]\n";
    }
    std::cout << "========================================================\n";

    return 0;
}

