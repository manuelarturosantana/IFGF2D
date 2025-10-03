#pragma once

#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <array>



#include "../Geometry/Patch.hpp"
#include "GreenFunctions.hpp"
#include "../IFGF_Source/BoxTree.h"

/*
Potential Future Improvements
     --> determine_patch_near_singular_points currently assumes a closed curve that does not
         nearly self intersect so that near singular points only come from adjecent patches
     --> I suspect most of the near singular points are very near the endpoints of the patches
         We may want to use Brent's method or perhaps just project to the patch enpoint to
         If the gss part becomes much to slow
    --> We may want two versions of the Green Function. One for sing/nearsing interactions that
        takes into account cancelation errors, and one for IFGF which doesn't need to account for such errors
    --> Ad hoc checking a point is actually on the boundayr or not
    --> *** Automatic selection of delta parameter
    --> *** Automatic determination of number of boxes
*/

/* 
Current assumptions to watch out for
    1. The number of points in each patch is the same. Will need to change the dermination
        of the near singular points if this changes. Actually will need each patch to store NP.
        Also will need to change IFGF checking if the sing/near sing points are in the neighborhood.
        Also init_sort_sing_point
    2. Currently the Morton box codes only go to 16 levels, which is about k = 512 we add a level
       everytime we double the wave number. The splitby2 and mergeFrom2 functions in IFGF_Source/Level.h 
       will need to be modified, as well as changing long longs to __uint128_t to make this work.
    3. Singularities only happen at the endpoints of curves.
    4. Open curves are always singular at both endpoints. If this is not the case, merge the parameterization
       into a single smooth parameterization.
    5. Corners are built from the junction of two open curves. ie there is not closed global parameterization
*/


/// @brief Class which computes the forward map of the greens function operator
/// @tparam Np Number of points to use per patch for integration
/// @tparam Formulation: Which Green's function formulation to use.
/// @tparam Nroot number of points to used in bounding box computation of patches
template <int Np, FormulationType Formulation, int Ps = 5, int Pang = 5, int Nroot = 10>
class ForwardMap {
    public:
        std::vector<Patch<Nroot>> patches_; 
        long long num_patches_;
        long long total_num_unknowns_;

        double delta_; // Rectangular Polar Near Singularity Parameter

        // Due to lifetime issues I am not going to store this as references
        // A vector of curve objects
        std::vector<std::unique_ptr<Curve>> curves_;
        // For each curve ii, store a vector of arrays where the first element of the array
        // says curve ii touches curve jj, and the second element tracks at which endpoint.
        // 0 for tlim1, or 1 for tlim2
        std::vector<std::vector<Junction>> curve_touch_tracker_;

        double p_; // Parameter in rectangular polar change of variables.
        double eta_ = 0.0; // Coupling parameter


        const int num_ns = 40; // Number of points to use for the near singular integration. Not a parameter to set for now
        
        std::vector<double> xs_; // x and y coordinates of all points on all patches.
        std::vector<double> ys_;
        std::vector<double> nxs_; // x and y coordinates of the unit normal vector
        std::vector<double> nys_;

        Eigen::ArrayXd fejer_nodes_;
        Eigen::ArrayXd fejer_weights_;

        Eigen::ArrayXd fejer_nodes_ns_;
        Eigen::ArrayXd fejer_weights_ns_;

        int GMRES_MAX_ITER_ = 200;
        double GMRES_TOLERANCE_ = 1e-6;

        // BoxTree<Ps, Pang> boxes_;

        // sort_sing_point_[i] returns a set tracking incides of points j, that are singular
        // or near singular to the point i. Here the indices i, j are in morton sorted order.
        std::vector<std::unordered_set<long long>> sort_sing_point_;
        

        
        ForwardMap(double delta, std::vector<std::unique_ptr<Curve>> curves, 
            const std::vector<std::vector<Junction>>& curve_touch_tracker, 
            double wavelengths_per_patch, 
            double patch_split_wavenumber, double p = 4);

        /// @brief Using a closed curve, computes the patches, and initializes their bounding boxes
        ///        Note, I know that you don't need to pass data members in, but passing in
        ///        curves and the tracker makes debugging easier.        
        /// @tparam N The number of points to use in the bounding box computation
        /// @param wavelengths_per_patch Number of wavelengths per patch
        /// @param wavenumber Wavenumber for the problem
        /// @return A vector of patches, which can be passed into the forward map
        void init_points_and_patches(
                double wavelengths_per_patch, double wavenumber);

        ///@brief Compute for each patch the point indices of the near singular points
        /// Note: To me this belongs here, and not in some geometry module because it depends
        /// on the number of discretization points in each patch
        void determine_patch_near_singular_points();
        
        /// Loop through the patches and compute the precomputations for each patch
        /// @param wavenumber The wavenumber for the solve. Not necessarily the same as for
        ///        init_points_and_patches, but should have real part <= the one in that function
        ///        to maintain accuracy.
        void compute_precomputations(std::complex<double> wavenumber);
        
        /// @brief Given a patch compute the precomputations corresponding to the weights.
        /// 
        /// @param patch The patch to integrate over
        /// @param tsing The in patch t-value corresponding to the singular point, or projected
        ///              point for the near singularity.
        /// @param xsing (also ysing) The x and y location of the singular (or near singular)
        ///        point. Necessary to pass in since for near singular point curve.x_t(t_sing) != xsing
        /// @param wave_number The wave number for this use of the solver.
        /// @param out_precomps  The vector or column of an eigenmatrix to be computed.
        /// @return The precomutations as a vector for each point
        Eigen::VectorXcd single_patch_point_compute_precomputations
            (const Patch<Nroot>& patch, double tsing, double xsing, double ysing, 
                 std::complex<double> wave_number, Eigen::Ref<Eigen::VectorXcd> out_precomps);

        /// @brief Multiply the density by the curve jacobian and integration weights
        ///        in preparation for using IFGF
        void compute_intensities(Eigen::VectorXcd& density);
        
        /// @brief Computes the non-singular interations in computing Ax without IFGF
        /// @param density The input density with intensities already computed
        /// @param wave_number The wave number of the problem
        /// @return A vector of the output
        Eigen::VectorXcd compute_non_singular_interactions_unacc(Eigen::VectorXcd& density, 
                std::complex<double> wave_number);
        
        /// @brief Returns a vector of the singular and near singular interactions given the density.
        Eigen::VectorXcd compute_sing_near_sing_interactions(Eigen::VectorXcd& density);

        /// @brief Compute the forward map Ax = b without IFGF
        Eigen::VectorXcd compute_Ax_unacc(Eigen::VectorXcd& density, std::complex<double> wave_number);
        void compute_Ax_unacc(const Eigen::VectorXcd& density, Eigen::VectorXcd& out, std::complex<double> wave_number);
        
        /// @brief Compute IFGF Precomputations and BoxTree set up
        /// @param wavenumber The wave number to solve at
        /// @param nlevels The number of levels to use in IFGF
        void precomps_and_setup(std::complex<double> wavenumber, int nlevels, BoxTree<Ps,Pang>& boxes); 
        
        /// @brief Initializes the tracking of singular/near singular points
        /// @param inverse The inverse sorting vector from Boxtree
        void init_sort_sing_point(const std::vector<long long>& inverse);


        /// @brief Compute the forward map Ax = b with IFGF 
        Eigen::VectorXcd compute_Ax_acc(Eigen::VectorXcd& density, BoxTree<Ps,Pang>& boxes);
        void compute_Ax_acc(const Eigen::VectorXcd& density, Eigen::VectorXcd& out, BoxTree<Ps,Pang>& boxes);
 
        Eigen::VectorXcd solve(const Eigen::VectorXcd& rhs, std::complex<double> wavenumber, 
            int nlevels=5,bool accelerated=true); 

        /// @brief Evaluate the layer potential via direct quadrature
        /// @param x_prop, y_prop x and y values to evaluate the layer potential at
        /// @param density The density to propigate. IMPORTANTLY passed as a copy so that
        ///                the weights can be multiplied before doing the propagation
        /// @param wavenumber The same wavenumber used to solve the equation
        Eigen::VectorXcd propagate_unacc(const std::vector<double> x_prop, 
            const std::vector<double> y_prop, Eigen::VectorXcd density, std::complex<double> wavenumber);

};


#include "ForwardMap.tpp"