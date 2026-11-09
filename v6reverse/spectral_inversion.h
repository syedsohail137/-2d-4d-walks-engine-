#ifndef SPECTRAL_INVERSION_H
#define SPECTRAL_INVERSION_H

#include <string>
#include <vector>

struct SpectralNode {
    int index;
    double riemann_zero_gamma;
    int prime;
    double phase_angle_rad;
    double critical_proximity;
    int target_determinant_d;
    std::string lattice_state;
    int a, b, c, d;
};

const std::vector<double> RIEMANN_GAMMAS = {
    14.13472514, 21.02203964, 25.01085758, 30.42487613, 32.93506159, 37.58617816,
    40.91871901, 43.32707328, 48.00515088, 49.77383248, 52.97032148, 56.44624770,
    59.34704400, 60.83177852, 65.11254405, 67.07981053, 69.54640171, 72.06715767,
    75.70469070, 77.14484007, 79.33737537, 82.91038085, 84.73549298, 87.42527461,
    89.24408101, 92.49162608, 95.38388277, 98.83119422, 101.44222071, 103.72553804
};

#endif // SPECTRAL_INVERSION_H
