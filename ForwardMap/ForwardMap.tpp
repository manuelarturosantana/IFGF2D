#pragma once
#include <chrono>
#include <iostream>
#include <array>

#include "ForwardMap.hpp"
#include "../utils/Quadrature.hpp"
#include "../utils/BrentsMethod.hpp"


template <int Np, FormulationType Formulation, int Nroot>
ForwardMap<Np, Formulation, Nroot>::ForwardMap(double delta, ClosedCurve& curve, double wavelengths_per_patch, 
    double patch_split_wavenumber, int near_singular_patch_est, double p_) :
        delta_(delta), curve_(curve), near_singular_patch_est_(near_singular_patch_est), p_(p)
{   
    // Precompute Quadrature nodes and weights for non-singular integration
    fejer_nodes_ = std::vector<double>(Np); fejer_weights_ = std::vector<double>(Np);
    fejerquadrature1(fejer_nodes_, fejer_weights_, Np);

    // Precopmute quadrature nodes and weights for singular integration
    fejer_nodes_ns_ = std::vector<double>(num_ns); fejer_weights_ns_ = std::vector<double>(num_ns);
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
    Eigen::ArrayXd cheb_left_nm1  =  Eigen::ArrayXd::Constant(num_ns, 1.0);  
    Eigen::ArrayXd cheb_right_nm1 =  Eigen::ArrayXd::Constant(num_ns, 1.0); 
    Eigen::ArrayXd cheb_left_n(num_ns);    Eigen::ArrayXd cheb_right_n(num_ns);
    Eigen::ArrayXd cheb_left_np1(num_ns);  Eigen::ArrayXd cheb_right_np1(num_ns);

    // Store the x values needed to update the Chebyshev recurrence
    Eigen::ArrayXd tl_cov(num_ns);  Eigen::ArrayXd tr_cov(num_ns);
    
    // Change of variables and curve jacobians. 
    Eigen::ArrayXd w_left(num_ns), w_right(num_ns), wp_left(num_ns), 
    wp_right(num_ns), cj_left(num_ns), cj_right(num_ns);;
    
    // Set up Green function evaluations on the left and right
    Eigen::VectorXcd  green_left(num_ns), green_right(np);

    const double a = patch.t1; const double b = patch.t2;
    const double xsing = curve_.xt(tsing); const double ysing = curve.yt(tsing);

    for (int ii = 0; ii < num_ns; ii++) {
        const double t_left  = (fejer_nodes_ns_(ii) + 1.0) / 2.0;
        const double t_right = (fejer_weights_ns_(ii) - 1.0) / 2.0;

        w_left(ii) =  w_cov(-t_left, p_);  wp_left(ii) = wp_cov(-t_left, p_);    
        w_right(ii) = w_cov(t_right, p_);  wp_right(ii) = wp_cov(t_right, p_);

        tl_cov(ii) = tsing - ((tsing - a) * w_left(ii));
        tr_cov(ii) = tsing + ((b - tsing) * w_right(ii));


        cheb_left_n(ii) = tl_cov(ii); cheb_right_n(ii) = tr_cov(ii);

        // Left Side
        const double xtl = curve.xt(ab2cd(tl_cov(ii) ,-1,1,a,b));
        const double ytl = curve.xt(ab2cd(tl_cov(ii) ,-1,1,a,b));

        Eigen::Vector2d normal_l = curve_.normal_t(tl_cov(ii) );
        cj_left(ii) = normal_l.norm();
        normal_l /= normal_l.norm();

        GF<Formulation>(xsing, ysing, xtl, ytl, normal_l(0), normal_l(1),
        eta_, wave_number, green_left.data() + ii);

        // Right Side
        const double xtr = curve.xt(ab2cd(tr_cov(ii),-1,1,a,b));
        const double ytr = curve.xt(ab2cd(tr_cov(ii),-1,1,a,b));

        Eigen::Vector2d normal_r = curve_.normal_t(tr_cov(ii));
        cj_right(ii) = normal_r.norm();
        normal_r /= normal_r.norm();

        GF<Formulation>(xsing, ysing, xtr, ytr, normal_r(0), normal_r(1),
        eta_, wave_number, green_right.data() + ii)
    }

    Eigen::ArrayXcd integrand_wo_cheb_l = ((tsing - a) / 2.0) * fejer_weights_ns_.array() * 
                            green_left * cj_left *  wp_left;
    Eigen::ArrayXcd integrand_wo_cheb_r = ((b - tsing) / 2.0) * fejer_weights_ns_.array() * 
                            green_right * cj_right * wp_right;

  

    // No need to element wise multiply by all ones.
    out_precomps(0) += integrand_wo_cheb_l.sum();
    out_precomps(0) += integrand_wo_cheb_r.sum();

    out_precomps(1) += (integrand_wo_cheb_l * cheb_left_n).sum();
    out_precomps(1) += (integrand_wo_cheb_r * cheb_right_n).sum();

    for (int ii = 2; ii < Np; ii++) {
        cheb_left_np1 = 2.0 * tl_cov  * cheb_left_n - cheb_left_nm1;
        cheb_right_np1 = 2.0 * tr_cov  * cheb_right_n - cheb_right_nm1;

        out_precomps(ii) += (integrand_wo_cheb_l * cheb_left_np1).sum();
        out_precomps(ii) += (integrand_wo_cheb_r * cheb_right_np1).sum();

        // Update the Chebyshev polynomials
        cheb_left_nm1.swap(cheb_left_n);   cheb_left_n.swap(cheb_left_np1);
        cheb_right_nm1.swap(cheb_right_n); cheb_right_n.swap(cheb_right_np1);
    }

}

// double near_sing_int(const std::function<double (double)>& f, double a, double b, double x_sing, 
//         int N, double p) {

//     Eigen::VectorXd nodes(N), weights(N); 

//     fejerquadrature1(nodes, weights, N);

//     Eigen::VectorXd w_left(N), w_right(N), wp_left(N), wp_right(N), fs_left(N), fs_right(N);
//     for (int ii = 0; ii < N; ii++) {
//         const double x_left  = (nodes(ii) + 1.0) / 2.0;
//         const double x_right = (nodes(ii) - 1.0) / 2.0;

//         w_left(ii) =  w_cov(-x_left, p);  wp_left(ii) = wp_cov(-x_left, p);    
//         w_right(ii) = w_cov(x_right, p);  wp_right(ii) = wp_cov(x_right,p);

//         fs_left(ii) = f(x_sing - ((x_sing - a) * w_left(ii)));
//         fs_right(ii) = f(x_sing + ((b - x_sing) * w_right(ii)));
//     }

   
//     fs_left = ((x_sing - a) / 2.0) * fs_left;
//     fs_right = ((b - x_sing) / 2.0) * fs_right;

//     return (weights.array() * fs_left.array() * wp_left.array()).sum() 
//         + (weights.array() * fs_right.array() * wp_right.array()).sum();

// }



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

