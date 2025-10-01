#include "../../Geometry/Curve.hpp"
#include <iostream>

int main() {
    // double wavenumber = 2;
    // double wavelengths_per_patch = 1.0;
    // Circle curve(0.0, 0.0, 1.0);

    double wavenumber = 8;
    double wavelengths_per_patch = 1.0;
    Line curve(1.0,1.0,10.0,12.0);
    // Kite curve;
    
    // Compute patch limits
    std::vector<double> patch_limits = curve.compute_patch_lims(wavelengths_per_patch, wavenumber);
    std::cout << "Patch limits size " << patch_limits.size() << std::endl;
    // Output the patch limits
    for (size_t i = 0; i < patch_limits.size(); ++i) {
        std::cout << "Patch " << i << ": [" << patch_limits[i] << ", ";
        if (i + 1 < patch_limits.size()) {
            std::cout << patch_limits[i + 1] << "]\n";
        } else {
            std::cout << "end]\n";
        }
    }

    // print for matlab readable format
    for (size_t i = 0; i < patch_limits.size(); ++i) {
        std::cout << patch_limits[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}