#include <iostream>
#include <functional>
#include <array>
#include <chrono>
#include<iomanip>

#include "../../complex_bessel-master/include/complex_bessel.h"

#include "../../ForwardMap/ForwardMap.hpp"
#include "../../Geometry/Curve.hpp"
#include "../../Geometry/Obstacles.hpp"
#include "../../utils/Rhs.hpp"
#include "../../utils/Chebyshev1d.hpp"

// File to test solution of the integral equation using GMRES
int main() {
    double delta = 0.1;
    double k     = 3;
    double wave_lengths_per_patch = 0.5;
    double x_0 = 0; double y_0 = 0; // Center of the point sources
    int nlevels = 6;
    // std::vector<double> x_test = {1.2}; std::vector<double> y_test = {1.5};
    std::vector<double> x_test = {0.17}; std::vector<double> y_test = {0.62};

    std::vector<std::unique_ptr<Curve>> curves;
    std::vector<std::vector<Junction>> curve_touch_tracker;

    // int nlevels = 8;

    // Circle test case
    // curves.emplace_back(std::make_unique<Circle>());
    // Square test case
    // std::vector<std::unique_ptr<Curve>> curves;
    // std::vector<std::vector<Junction>> curve_touch_tracker;
    // make_rect(curves, curve_touch_tracker);

    // Straight Line test case
    curves.emplace_back(std::make_unique<Line>(-1, -0.2,1,-0.2));


    constexpr int num_points = 20;

    ForwardMap<num_points, FormulationType::SingleLayer, 5, 5> FM(delta, std::move(curves), 
        curve_touch_tracker, wave_lengths_per_patch, k, 6);
    std::cout << "Num unknowns " << FM.total_num_unknowns_ << std::endl;

    // Eigen::VectorXcd RHS = point_source(FM.xs_, FM.ys_, x_0, y_0, k);
    Eigen::VectorXcd RHS = test_case_straight_line(FM.xs_, FM.ys_);

    // Eigen::VectorXcd dens = FM.solve(RHS, k, nlevels, false);
    Eigen::VectorXcd dens = FM.solve(RHS, k, nlevels, true);

    Eigen::VectorXcd sol_app = FM.propagate_unacc(x_test, y_test, dens, k);

    double distance = std::sqrt(
        std::pow(x_test[0] - x_0, 2) + std::pow(y_test[0] - y_0, 2)
    );


    // std::complex<double> sol_true = sp_bessel::hankelH1(0, k * distance);

    std::complex<double> sol_true(0.02788626934981090, - 0.75932847390327920);

    std::cout << "Error in Solution " << std::abs(sol_true - sol_app(0)) << std::endl;
    
    return 0;

}