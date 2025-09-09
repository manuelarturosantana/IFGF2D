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
    double singular_tval = 0;
    // double singular_tval = 0.392699;

    // Circle curve;
    Kite curve;
    constexpr int num_points = 10;
    ForwardMap<num_points, FormulationType::SingleLayer, num_points> FM(delta, curve, wave_lengths_per_patch, k);

        // Prints the matlab code to plot the patches and bounding boxes
      for (size_t i = 0; i < 1; i++) {
        std::cout 
        << "patch = curve.x_t(linspace( " << FM.patches_[i].t1 << "," << FM.patches_[i].t2 <<")); \n";            
    }
    std::cout << std::endl;

    // TODO: This appears to work, but only for 5-6 digits matching Matlab's function
    // I believe this is due to cancellation errors in the Sing Layer Potential
    // We may fix this in the future, but since we are using IFGF 5-6 digits may suffice
    Eigen::IOFormat FullPrecisionFmt(Eigen::FullPrecision);
    Eigen::VectorXcd out(num_points);
    double xsing = curve.xt(singular_tval); double ysing = curve.yt(singular_tval);
    FM.single_patch_point_compute_precomputations(FM.patches_[0], singular_tval, xsing, ysing, k, out);
    std::cout << "Precomputations \n" << out.format(FullPrecisionFmt) << std::endl;
    


    return 0;
}