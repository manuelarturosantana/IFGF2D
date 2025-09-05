#pragma once

#include <vector>
#include <memory>


#include "../Geometry/Patch.hpp"
#include "GreenFunctions.hpp"

/*
Potential Future Improvements
     --> determine_patch_near_singular_points currently assumes a closed curve that does not
         nearly self intersect so that near singular points only come from adjecent patches
     --> I suspect most of the near singular points are very near the endpoints of the patches
         We may want to use Brent's method or perhaps just project to the patch enpoint to
         If the gss part becomes much to slow
     --> Rather than GSS for now I am just going to do the change of variables at the endpoints
         of the patch. As long as it is not nearly self intersecting this should work.

*/

/* 
Current assumptions to watch out for
    1. The number of points in each patch is the same. Will need to change the dermination
        of the near singular points if this changes

*/



/// @brief Class which computes the forward map of the greens function operator
/// @tparam Np Number of points to use per patch for integration
/// @tparam Formulation: Which Green's function formulation to use.
/// @tparam Nroot number of points to used in bounding box computation of patches
template <int Np, FormulationType Formulation, int Nroot = 10>
class ForwardMap {
    public:
        std::vector<Patch<Nroot>> patches_; 
        long long num_patches_;
        long long total_num_unknowns_;

        ClosedCurve& curve_;

        double delta_; // Rectangular Polar Near Singularity Parameter
        int near_singular_patch_est_; // Assume all near singular points are in patches only 1 away from this
        double p_; // Parameter in rectangular polar change of variables.
        double gss_tol_;

        const num_ns_integral_points = 20; // Number of points to use on each side of the near 
                                           // singular integration. Not a parameter to set for now
        
        std::vector<double> xs_; // x and y coordinates of all points on all patches.
        std::vector<double> ys_;

        std::vector<double> fejer_nodes_;
        std::vector<double> fejer_weights_;

        
        ForwardMap(double delta, ClosedCurve& curve, double wavelengths_per_patch, 
            double patch_split_wavenumber, int near_singular_patch_est_ = 1, double p_ = 6);

        /// @brief Using a closed curve, computes the patches, and initializes their bounding boxes
        /// @tparam N The number of points to use in the bounding box computation
        /// @param curve The curve to be used
        /// @param wavelengths_per_patch Number of wavelengths per patch
        /// @param wavenumber Wavenumber for the problem
        /// @return A vector of patches, which can be passed into the forward map
        void init_points_and_patches(ClosedCurve& curve, 
                double wavelengths_per_patch, double wavenumber);

        ///@brief Compute for each patch the point indices of the near singular points
        /// Note: To me this belongs here, and not in some geometry module because it depends
        /// on the number of discretization points in each patch
        void determine_patch_near_singular_points();
        

        

        // Loop through the patches and compute the precomputations for each patch
        void compute_precomputations();
        
        /// @brief Given a patch compute the precomputations corresponding to the weights.
        /// Works if the point is not on the endpoints 
        /// @param patch The patch to integrate over
        /// @param tval The in patch t-value of the integration point
        /// @return An eigen vector corresponding to the weights.
        Eigen::VectorXcd single_patch_point_mid_compute_precomputations(const Patch<Np>& patch, double tsing);


        // Other functions to implement:
        // Compute intensities: Take in the density and multiply by the correct weights
        // compute_forward_map: Send x -> Ax by (Take in k as an argument)
            // Computing the Singular / Near singular interactions
            // Computing the intensities and doing the direct sum
            // Eventually going to include IFGF
};


#include "ForwardMap.tpp"