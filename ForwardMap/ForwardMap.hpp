#pragma once

#include <vector>
#include <memory>


#include "../Geometry/Patch.hpp";

/// @brief Class which computes the forward map of the greens function operator
/// @tparam Np Number of points to use per patch for integration
/// @tparam Nroot number of points to used in bounding box computation of patches
template <int Np, int Nroot = 10>
class ForwardMap {
    public:
        std::vector<Patch<Nroot>> patches_; 
        long long num_patches_;
        long long total_num_unknowns;

        double delta_; // Rectangular Polar Near Singularity Parameter
        
        std::vector<double> xs_; // x and y coordinates of all points on all patches.
        std::vector<double> ys_;

        std::vector<double> fejer_nodes;
        std::vector<double> fejer_weights;

        
        ForwardMap(double delta, ClosedCurve& curve, double wavelengths_per_patch, double patch_split_wavenumber);

        /// @brief Using a closed curve, computes the patches, and initializes their bounding boxes
        /// @tparam N The number of points to use in the bounding box computation
        /// @param curve The curve to be used
        /// @param wavelengths_per_patch Number of wavelengths per patch
        /// @param wavenumber Wavenumber for the problem
        /// @return A vector of patches, which can be passed into the forward map
        std::vector<Patch<Nroot>> init_points_and_patches(ClosedCurve& curve, 
                double wavelengths_per_patch, double wavenumber);


        // Because the xs, and ys are dependent on the number of points per patch,
        // this function belongs in the Forward map
        void compute_near_singular_points();

        // Loop through the patches and compute the precomputations for each patch
        void compute_precomputations();

        Eigen::MatrixXcd single_patch_compute_precomputations(const Patch<N>& patch, int npoints, double p);

        // Other functions to implement:
        // Compute intensities: Take in the density and multiply by the correct weights
        // compute_forward_map: Send x -> Ax by (Take in k as an argument)
            // Computing the Singular / Near singular interactions
            // Computing the intensities and doing the direct sum
            // Eventually going to include IFGF
};



// template <int N = 10>
// std::vector<Patch<N>> ForwardMap::initialze_patches(ClosedCurve& curve, double wavelengths_per_patch, double wavenumber) {
//     std::vector<double> patch_lims = curve.compute_patch_lims(wavelengths_per_patch, wavenumber);
//     std::vector<Patch<N>> patches; patches.reserve(patch_lims.size() - 1);

//     for (size_t ii = 0; ii < patch_lims.size() - 1; ii++) {
//         patches.emplace_back(patch_lims[ii], patch_lims[ii + 1], curve, 0.1);
//         patches[i].comp_bounding_box(); 
//     }
// }

#include "ForwardMap.tpp"