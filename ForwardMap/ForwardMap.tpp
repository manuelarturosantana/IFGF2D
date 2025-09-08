#pragma once
#include <chrono>
#include <iostream>
#include <array>

#include "ForwardMap.hpp"
#include "../utils/Quadrature.hpp"
#include "../utils/BrentsMethod.hpp"


template <int Np, FormulationType Formulation, int Nroot>
ForwardMap<Np, Formulation, Nroot>::ForwardMap(double delta, ClosedCurve& curve, double wavelengths_per_patch, 
    double patch_split_wavenumber, int near_singular_patch_est, double p) :
        delta_(delta), curve_(curve), near_singular_patch_est_(near_singular_patch_est), p_(p)
{   
    // Precompute Quadrature nodes and weights for non-singular integration
    fejer_nodes_ = std::vector<double>(Np); fejer_weights_ = std::vector<double>(Np);
    fejerquadrature1(fejer_nodes_, fejer_weights_, Np);

    // Precopmute quadrature nodes and weights for singular integration
    fejer_nodes_ns_ = Eigen::ArrayXd(num_ns); fejer_weights_ns_ =  Eigen::ArrayXd(num_ns);
    fejerquadrature1(fejer_nodes_ns_, fejer_weights_ns_, num_ns);

    
    // Compute all points used in the discretization and bounding boxes.
    init_points_and_patches(curve, wavelengths_per_patch, patch_split_wavenumber);

    auto start = std::chrono::high_resolution_clock::now();
    determine_patch_near_singular_points();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to determine near singular points: " << elapsed.count() << std::endl;
   

}

template <int Np, FormulationType Formulation, int Nroot>
void ForwardMap<Np, Formulation, Nroot>::init_points_and_patches(ClosedCurve& curve, 
    double wavelengths_per_patch, double wavenumber) {

    std::vector<double> patch_lims = curve.compute_patch_lims(wavelengths_per_patch, wavenumber);
    
    num_patches_ = patch_lims.size() - 1;
    total_num_unknowns_ = num_patches_ * Np;

    patches_.clear(); patches_.reserve(num_patches_);

    xs_.clear(); ys_.clear(); 
    xs_.reserve(num_patches_ * total_num_unknowns_); 
    ys_.reserve(num_patches_ * total_num_unknowns_);

    for (size_t ii = 0; ii < num_patches_; ii++) {
        patches_.emplace_back(patch_lims[ii], patch_lims[ii + 1], curve, delta_);
        patches_[ii].comp_bounding_box();

        // Put points in a vector on the curve
        for (int jj = 0; jj < Np; jj++) { 
            double t_curr = ab2cd(fejer_nodes_[jj], -1, 1, patch_lims[ii], patch_lims[ii + 1]);
            xs_.push_back(curve.xt(t_curr));
            ys_.push_back(curve.yt(t_curr));
        }
    }
}

template <int Np, FormulationType Formulation, int Nroot>
void ForwardMap<Np, Formulation, Nroot>::determine_patch_near_singular_points() {
    for (int ii = 0; ii < num_patches_; ii++) {
    // Pointer to the current patch to test the near singular points.
    Patch<Nroot>* p_tp = &patches_[ii];
        for (int jj = ii - near_singular_patch_est_; jj <= ii + near_singular_patch_est_; jj++) {

            int test_patch_ind;

            // Look at neighboring patches, including endpoints.
            if (jj == ii) {
                continue;
            } else if (jj < 0) {
                test_patch_ind = num_patches_ + jj;
                // >= accounts for zero indexing
            } else if (jj >= num_patches_) {
                test_patch_ind = jj - num_patches_;
            } else {
                test_patch_ind = jj;
            }

            // Loop over the points in the patch, and check if it is near singular
            // We use the global and not patchwise ordering so we can store which points
            // are near singular to a certain patch
            int tp_val_start = test_patch_ind * Np;

            for (int kk = tp_val_start; kk < tp_val_start + Np; kk++) {
                // DEBUG
                // std::cout << "Bounding box computation for " << kk << std::endl;

                // See if it is in the bounding box
                double xp = xs_[kk];
                double yp = ys_[kk];
               
                if (p_tp->bounding_box_.is_inside(xp, yp)) {
                    
                    // Lambda returning the distance of f to a point on the curve
                    auto f = [xp, yp, this] (double t) -> double {
                        double xt = this->curve_.xt(t);
                        double yt = this->curve_.yt(t);
                        return std::sqrt(std::abs( 
                            (xt - xp) * (xt - xp) + (yt - yp) * (yt - yp)
                           ));
                    };

                    if (f(p_tp->t1) < delta_) {

                        p_tp->near_singular_point_indices_.push_back(kk);
                        p_tp->near_singular_point_ts_.push_back(p_tp->t1);

                    } else if (f(p_tp->t2) < delta_) {

                        p_tp->near_singular_point_indices_.push_back(kk);
                        p_tp->near_singular_point_ts_.push_back(p_tp->t2);

                    } else {

                        std::array<double,2> out = brent_find_minima(f, p_tp->t1, p_tp->t2);
             
                        // Check if the point in the bounding box is near singular.
                        if (out[1] < delta_) {
                            // // DEBUG
                            // std::cout << " Found a near singular point via brent_method" << std::endl;
                            // std::cout << out[0] << "  " << out[1] << std::endl;
                            // std::cout << "Patch lims" <<  patches_[test_patch_ind].t1 << " " << patches_[test_patch_ind].t2;
                            p_tp->near_singular_point_indices_.push_back(kk);
                            p_tp->near_singular_point_ts_.push_back(out[0]);
                        }

                    }
                    
                }

            }

        }

    }
}

template <int Np, FormulationType Formulation, int Nroot>
Eigen::VectorXcd ForwardMap<Np, Formulation, Nroot>::single_patch_point_mid_compute_precomputations
    (const Patch<Np>& patch, double tsing, double wave_number) {

    // Setup output vector
    Eigen::VectorXcd out_precomps = Eigen::VectorXcd::Zero(Np);

    // Set up vector of Cheb left and Cheb right
    Eigen::ArrayXd cheb_nm1  =  Eigen::ArrayXd::Constant(num_ns, 1.0);  
    Eigen::ArrayXd cheb_n(num_ns);    Eigen::ArrayXd cheb_right_n(num_ns);
    Eigen::ArrayXd cheb_np1(num_ns);  Eigen::ArrayXd cheb_right_np1(num_ns);

    // Store the change of variables evaluated at the cheb points
    Eigen::ArrayXd xi_alpha(num_ns);
    
    // Change of variablesfrom [t1,t2] to [-1,1] cov and curve jacobians. 
    Eigen::ArrayXd g_cov(num_ns), curve_jac(num_ns);
    
    // Set up Green function evaluations on the left and right
    Eigen::ArrayXcd  green(num_ns);
    
    // Cartesean coordinates of the singular point.
    const double xsing = curve_.xt(tsing); const double ysing = curve_.yt(tsing);
    
    // Location in [-1,1] of the singularity
    double alpha = ab2cd(tsing,-1.0, 1.0, patch.t1, patch.t2);

    // 0 for middle. -1 for left, 1 for right.
    int sing_loc;
    
    // This check to see if it is on the enpoint is ad hoc for now. Perhaps it needs to be tighter
    if (std::abs(-1.0 - alpha) < 1e-5) { 
        sing_loc = -1;
    } else if (std::abs(1.0 - alpha) < 1e-5) {
        sing_loc = 1;
    } else {
        sing_loc = 0;
    }
    
    // Precompute all values which don't depend on the Chebyshev polynomials.
    for (int ii = 0; ii < num_ns; ii++) {
        
        // Compute the singularity resolving change of variables.
        if (sing_loc == 0) {
            xi_alpha(ii) = xi_cov_center(fejer_nodes_ns_(ii), alpha, p_);
        } else if (sing_loc == -1) {
            xi_alpha(ii) = xi_cov_left(fejer_nodes_ns_(ii), p_);
        } else {
            xi_alpha(ii) = xi_cov_right(fejer_nodes_ns_(ii), p_);
        }

        cheb_left_n(ii) = xi_alpha(ii);

        // Compute the change of variables from [-1,1] to t1, t2;
        double g_xi_alpha = ab2cd(xi_alpha(ii), -1.0, 1.0, patch.t1, patch.t2);

        const double xt = curve_.xt(g_xi_alpha);
        const double yt = curve_.xt(g_xi_alpha);

        // Compute the curve Jacobian
        Eigen::Vector2d normal = curve_.normal_t(g_xi_alpha);
        curve_jac(ii) = normal.norm();

        // Compute the Green's Function
        normal       /= curve_jac(ii);
        green(ii) = GF<Formulation>(xsing, ysing, xt, yt, normal(0), normal(1),
        eta_, wave_number);
    }

    Eigen::ArrayXcd integrand_wo_cheb = 
        ((patch.t2 - patch.t1) / (2.0)) * green * curve_jac * xi_alpha * fejer_weights_ns_;

    // No need to element wise multiply by all ones.
    out_precomps(0) = integrand_wo_cheb.sum();

    out_precomps(1) = (integrand_wo_cheb * cheb_n).sum();

    for (int ii = 2; ii < Np; ii++) {
        cheb_np1 = 2.0 * xi_alpha  * cheb_n - cheb_nm1;

        out_precomps(ii) += (integrand_wo_cheb * cheb_np1).sum();

        // Update the Chebyshev polynomials
        cheb_left_nm1.swap(cheb_n);   
        cheb_left_n.swap(cheb_np1);
    }

    return out_precomps;

}

template <int Np, FormulationType Formulation, int Nroot>
void ForwardMap<Np, Formulation, Nroot>::compute_precomputations() {
    // Loop over patches
        // Init Patch eigen matrix of size npoints_rows x n_points + n_near singualr poitns columsn
        // Loop over x values 
            // Set the colum with an function which computes the precomputations
            // matrix.col = ()...
        // Loop over near singular points
            // If distance of projection point to endpoint is < 1e-3
            // Loop over Cheb polynomials
            // Do a once sided change of variables
            // Loop over Cheb polynomials
            // otherwise do the two sided change of variables

}

