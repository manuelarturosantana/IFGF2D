#include <iostream> 
// for timing
#include <vector>
#include <chrono>
#include <memory>

#include "../../utils/Chebyshev1d.hpp"
#include "../../Geometry/ClosedCurve.hpp"
#include "../../Geometry/Patch.hpp"

int main() {
    double k = 100.0;
    Kite kite;

    kite.set_integration_points();
    std::vector<double> patch_lims = kite.compute_patch_lims(1.0, k);
    std::vector<Patch<10>> patches;
    kite.integration_data_cleanup();

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < patch_lims.size() - 1; i++) {
        // std::make_unique creates a unique kite on the heap
        // patches.emplace_back(patch_lims[i], patch_lims[i + 1], std::make_unique<Kite>(), 0.1);
        patches.emplace_back(patch_lims[i], patch_lims[i + 1], kite, 0.1);
        patches[i].comp_bounding_box();    
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to compute bounding boxes for " << patches.size() 
              << " patches: " << elapsed.count() << " seconds" << std::endl;

    // Prints the matlab code to plot the patches and bounding boxes
    //   for (size_t i = 0; i < patch_lims.size() - 1; i++) {
    //     std::cout 
    //     << "patch = curve.x_t(linspace( " << patches[i].t1 << "," << patches[i].t2 <<")); \n" 
    //     << "plot(patch(1,:), patch(2,:)) \n" 
    //         << "x = ["  << patches[i].bounding_box_.Ax << ", " << patches[i].bounding_box_.Bx  
    //              << ", " << patches[i].bounding_box_.Cx << ", " << patches[i].bounding_box_.Dx << ", " 
    //              << patches[i].bounding_box_.Ax << "];\n"
    //         << "y = ["  << patches[i].bounding_box_.Ay << ", " << patches[i].bounding_box_.By  
    //              << ", " << patches[i].bounding_box_.Cy << ", " << patches[i].bounding_box_.Dy << ", "
    //              << patches[i].bounding_box_.Ay<< "];\n" 
    //         << "plot(x, y) \n";
            
    // }
    // std::cout << std::endl;

    return 0;
}


// // Some simple root finding tests
// // auto f = [](double x) { return std::cos(3.0 * x) + x; };  
// //double a = 0; double b = 1;
// auto f = [](double x) { return std::sin(x); };
// double a = 0; double b = 2 * M_PI;
// // Test case from approximation theory and pratice book.
// // auto f = [](double x) { return x * (x - (1.0/4.0)) * (x - 1.0 / 2.0); }; 
// // double a = -1; double b = 1;

// auto start = std::chrono::high_resolution_clock::now();
// std::vector<double> roots = Cheb1D::cheb_root_finder<10>(f, a, b);
// auto end = std::chrono::high_resolution_clock::now();
// std::chrono::duration<double> elapsed = end - start;
// std::cout << "Time taken for root finding: " << elapsed.count() << " seconds" << std::endl;

// // Print the roots
// std::cout << "Roots found in [" << a << ", " << b << "]:" << std::endl;
// for (const auto& root : roots) {
//     std::cout << root << std::endl;
// }