#include <iostream>
#include <functional>
#include <array>

#include "../../ForwardMap/ForwardMap.hpp"
#include "../../Geometry/ClosedCurve.hpp"
#include "../../utils/BrentsMethod.hpp"

// Appears to be working :).
// Can be used to test individual patches.
int main() {
    double delta = 0.1;
    double k     = 5;
    double wave_lengths_per_patch = 1.0;

    // Brents method test
    // auto f = [](double x) {
    //     return (x + 3) * (x - 1) * (x - 1);
    // };

    // double tlim1 = -1.0;
    // double tlim2 = 3;

    // std::array<double, 2> out = brent_find_minima(f, -1.0, 3.0);
    // std::cout << "Min = " << out[0] << " f(min) = " << out[1] << std::endl;

    
    Kite kite;
    constexpr int num_points = 10;
    ForwardMap<num_points, FormulationType::SingleLayer, num_points> FM(delta, kite, wave_lengths_per_patch, k);

    // Prints the matlab code to plot the patches and bounding boxes
    //    for (size_t i = 0; i < 1; i++) {
      for (size_t i = 0; i < FM.num_patches_; i++) {
        std::cout 
        << "patch = curve.x_t(linspace( " << FM.patches_[i].t1 << "," << FM.patches_[i].t2 <<")); \n" 
        << "plot(patch(1,:), patch(2,:)) \n" 
            << "x = ["  << FM.patches_[i].bounding_box_.Ax << ", " << FM.patches_[i].bounding_box_.Bx  
                 << ", " << FM.patches_[i].bounding_box_.Cx << ", " << FM.patches_[i].bounding_box_.Dx << ", " 
                 << FM.patches_[i].bounding_box_.Ax << "];\n"
            << "y = ["  << FM.patches_[i].bounding_box_.Ay << ", " << FM.patches_[i].bounding_box_.By  
                 << ", " << FM.patches_[i].bounding_box_.Cy << ", " << FM.patches_[i].bounding_box_.Dy << ", "
                 << FM.patches_[i].bounding_box_.Ay<< "];\n" 
            << "plot(x, y) \n";
            
    }
    std::cout << std::endl;

    std::cout << "x = [";
    for (int i = 0; i < FM.total_num_unknowns_; i++) {
        std::cout << FM.xs_[i] << "\n";
    }
    std::cout << "];" << std::endl;

    std::cout << "y = [";
    for (int i = 0; i < FM.total_num_unknowns_; i++) {
        std::cout << FM.ys_[i] << "\n";
    }
    std::cout << "];" << std::endl;
    std::cout << "plot(x,y,'*')" << std::endl;


    // for (int pind = FM.num_patches_  - 1 ; pind < FM.num_patches_; pind++) {
    for (int pind = 0; pind < FM.num_patches_; pind++) {
        for (size_t ns_ind = 0; ns_ind < FM.patches_[pind].near_singular_point_indices_.size(); ns_ind++) {
            double tval = FM.patches_[pind].near_singular_point_ts_[ns_ind];
 
            std::cout << 
            " tval = " << tval << "\n"
            << "xp = curve.x_t(tval)" << "\n"
            << "plot(xp(1), xp(2),'^')" << std::endl;
        }

    }




    return 0;
}