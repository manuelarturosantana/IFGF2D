#include <iostream>
#include <functional>
#include <array>
#include <chrono>
#include<iomanip>

#include "../../complex_bessel-master/include/complex_bessel.h"

#include "../../ForwardMap/ForwardMap.hpp"
#include "../../Geometry/ClosedCurve.hpp"
#include "../../utils/Rhs.hpp"
#include "../../utils/Chebyshev1d.hpp"

// File to test solution of the integral equation using GMRES
int main() {
    double delta = 0.01;
    double k     = M_PI * 10;
    double wave_lengths_per_patch = 1;
    double x_0 = 0.1; double y_0 = 0.1; // Center of the point sources

    std::vector<double> x_test = {1.2}; std::vector<double> y_test = {1.5};

    // int nlevels = 8;

    Circle circle;
    constexpr int num_points = 10;

    ForwardMap<num_points, FormulationType::SingleLayer, 5, 5> FM(delta, circle, wave_lengths_per_patch, k);
    std::cout << "Num unknowns " << FM.total_num_unknowns_ << std::endl;


    Eigen::VectorXcd RHS = point_source(FM.xs_, FM.ys_, k, x_0, y_0);

    Eigen::VectorXcd dens = FM.solve_unacc(RHS, k);

    Eigen::VectorXcd sol_app = FM.propagate_unacc(x_test, y_test, dens, k);

    double distance = std::sqrt(
        std::pow(x_test[0] - x_0, 2) + std::pow(y_test[0] - y_0,2)
    );

    std::complex<double> c_unit(0, 1.0);

    std::complex<double> sol_true = (c_unit / 4.0) * sp_bessel::hankelH1(0, k * distance);

    std::cout << "Error in Solution " << std::abs(sol_true - sol_app(0)) << std::endl;
    
    return 0;

}