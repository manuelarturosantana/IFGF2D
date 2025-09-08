#include <iostream>

#include <iostream>
#include <functional>
#include <array>

#include <Eigen/Dense>

#include "../../ForwardMap/ForwardMap.hpp"
#include "../../Geometry/ClosedCurve.hpp"


int main() {

    double delta = 0.1;
    double k     = 5;
    double wave_lengths_per_patch = 1.0;
    double singular_tval = 0.2;

    Circle circle;
    // Kite kite;
    constexpr int num_points = 10;
    ForwardMap<num_points, FormulationType::SingleLayer, num_points> FM(delta, circle, wave_lengths_per_patch, k);

        // Prints the matlab code to plot the patches and bounding boxes
      for (size_t i = 0; i < 1; i++) {
        std::cout 
        << "patch = curve.x_t(linspace( " << FM.patches_[i].t1 << "," << FM.patches_[i].t2 <<")); \n";            
    }
    std::cout << std::endl;

    // TODO: This appears to work, but only for 5-6 digits matching Matlab's function
    // Could be a matlab thing, or something with my function. This could be an issue when testing 
    // convergence. Something to check if the overall method doesn't end up as accurate as we hope
    Eigen::IOFormat FullPrecisionFmt(Eigen::FullPrecision);
    Eigen::VectorXcd out = FM.single_patch_point_compute_precomputations(FM.patches_[0], singular_tval, k);
    std::cout << "Precomputations \n" << out.format(FullPrecisionFmt) << std::endl;
    


    return 0;
}