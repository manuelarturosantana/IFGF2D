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


// Script to Test the accuracy of the forward map on an eigenfunction test.
//  Eigen::IOFormat customFormat(Eigen::FullPrecision, 0, ", ", "\n", "[", "]", "[", "]");
// std::cout << vec.format(customFormat) << std::endl;
int main() {
    double delta = 0.01;
    double k     = M_PI * 300;
    double wave_lengths_per_patch = 1;
    int m = 3;
    int nlevels = 8;

    Circle circle;
    constexpr int num_points = 10;

    ForwardMap<num_points, FormulationType::SingleLayer, 5, 5> FM(delta, circle, wave_lengths_per_patch, k);

    std::cout << "Total number of points: " << FM.total_num_unknowns_ << std::endl;

    Eigen::VectorXcd RHS = circle_eigenfunction(FM.xs_, FM.ys_, m);
 
    std::cout << "Num unknowns " << FM.total_num_unknowns_ << std::endl;
    std::cout << "Num Patches " << FM.num_patches_ << std::endl;
    std::cout << "Size of a patch " << FM.patches_[0].point_t_vals_.size() << std::endl;


    auto start = std::chrono::high_resolution_clock::now();
    BoxTree<5,5> boxes;
    FM.precomps_and_setup(k, nlevels, boxes);
    auto end   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to run precomps and setup: " << elapsed.count() << std::endl;

    
    start = std::chrono::high_resolution_clock::now();
    Eigen::VectorXcd Ax_unacc = FM.compute_Ax_unacc(RHS,k);
    end   = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time taken to compute Ax Unacc " << elapsed.count() << std::endl;

    RHS = circle_eigenfunction(FM.xs_, FM.ys_, m);

    start = std::chrono::high_resolution_clock::now();
    Eigen::VectorXcd Ax_acc = FM.compute_Ax_acc(RHS,boxes);
    end   = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time taken to compute Ax Acc " << elapsed.count() << std::endl;


    std::complex<double> c_unit(0.0, 1.0);
    std::complex<double> eval = ((c_unit * M_PI) / 2.0) * sp_bessel::besselJ(m,k) 
                * sp_bessel::hankelH1(m, k);

    std::cout << "Eigenvalue: " << eval << std::endl;
    // The RHS gets changed :(
    RHS = circle_eigenfunction(FM.xs_, FM.ys_, m);

    Eigen::VectorXcd diff = Ax_unacc - (eval * RHS);
    std::cout << "Error in the Eigenfunction Test Un-acc: " << diff.norm() << std::endl;
   
    diff = Ax_acc - (eval * RHS);
    std::cout << "Error in the Eigenfunction Test acc: " << diff.norm() << std::endl;

    return 0;

}