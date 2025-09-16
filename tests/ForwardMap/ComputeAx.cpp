#include <iostream>
#include <functional>
#include <array>
#include <chrono>

#include "../../complex_bessel-master/include/complex_bessel.h"

#include "../../ForwardMap/ForwardMap.hpp"
#include "../../Geometry/ClosedCurve.hpp"
#include "../../utils/Rhs.hpp"

// Script to Test the accuracy of the forward map on an eigenfunction test.

int main() {
    double delta = 0.1;
    double k     = 2;
    double wave_lengths_per_patch = 0.5;
    int m = 3;

    Circle circle;
    constexpr int num_points = 20;
    ForwardMap<num_points, FormulationType::SingleLayer, num_points> FM(delta, circle, wave_lengths_per_patch, k);

    auto start = std::chrono::high_resolution_clock::now();
    FM.compute_precomputations(std::complex<double>(k));
    auto end   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to compute precomputations: " << elapsed.count() << std::endl;

    Eigen::VectorXcd RHS = circle_eigenfunction(FM.xs_, FM.ys_, m);
    //DEBUG 
    std::cout << "Num unknowns " << FM.total_num_unknowns_ << std::endl;
    std::cout << "Num Patches " << FM.num_patches_ << std::endl;
    std::cout << "Size of a patch " << FM.patches_[0].point_t_vals_.size() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    Eigen::VectorXcd Ax = FM.compute_Ax_unacc(RHS,k);
    end   = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time taken to compute Ax " << elapsed.count() << std::endl;

    int test_patch = 0;
    std::cout << "Patch lims " << FM.patches_[test_patch].t1 << " " << FM.patches_[test_patch].t2 << std::endl;

    std::complex<double> c_unit(0.0, 0.1);
    std::complex<double> eval = ((c_unit * M_PI) / 2.0) * sp_bessel::besselJ(m,k) 
                * sp_bessel::hankelH1(m, k);
    
    Eigen::VectorXcd diff = Ax - (eval * RHS);

    std::cout << "Error in the Eigenfunction Test: " << diff.norm() << std::endl;

    return 0;

}