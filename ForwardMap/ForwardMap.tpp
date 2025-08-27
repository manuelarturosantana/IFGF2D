#pragma once
#include <chrono>
#include <iostream>
#include <array>

#include "ForwardMap.hpp"
#include "../utils/Quadrature.hpp"
#include "../utils/BrentsMethod.hpp"


template <int Np, int Nroot>
ForwardMap<Np, Nroot>::ForwardMap(double delta, ClosedCurve& curve, double wavelengths_per_patch, 
    double patch_split_wavenumber, int near_singular_patch_est) :
        delta_(delta), curve_(curve), near_singular_patch_est_(near_singular_patch_est)
{   
    // Precompute Quadrature nodes and weights
    fejer_nodes_ = std::vector<double>(Np); fejer_weights_ = std::vector<double>(Np);
    fejerquadrature1(fejer_nodes_, fejer_weights_, Np);
    
    // Compute all points used in the discretization and bounding boxes.
    init_points_and_patches(curve, wavelengths_per_patch, patch_split_wavenumber);

    auto start = std::chrono::high_resolution_clock::now();
    determine_patch_near_singular_points();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to determine near singular points: " << elapsed.count() << std::endl;
   

}

template <int Np, int Nroot>
void ForwardMap<Np, Nroot>::init_points_and_patches(ClosedCurve& curve, 
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

template <int Np, int Nroot>
void ForwardMap<Np, Nroot>::determine_patch_near_singular_points() {
    for (int ii = 0; ii < num_patches_; ii++) {
        for (int jj = ii - near_singular_patch_est_; jj < ii + near_singular_patch_est_; jj++) {

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
            int tp_val_start = test_patch_ind * Np;

            // Pointer to the test patch
            Patch<Nroot>* p_tp = &patches_[test_patch_ind];

            for (int kk = tp_val_start; kk < tp_val_start + Np; kk++) {
                // DEBUG
                // std::cout << "Bounding box computation for " << kk << std::endl;

                // See if it is in the bounding box
                double xp = xs_[kk];
                double yp = ys_[kk];
                if (p_tp->bounding_box_.is_inside(xp, yp)) {

                    // TO FIX:
                    //  Rather than GSS test distance by just measuring
                    //  The distance to the endpoint of an interval...
                    
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

                        // std::cout << "attempting to find minima via brents method" << std::endl;
                        std::array<double,2> out = brent_find_minima(f, p_tp->t1, p_tp->t2);

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

