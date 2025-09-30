#pragma once
#include <chrono>
#include <iostream>
#include <array>
#include <iomanip>

#include "opGMRES.hpp"

#include "ForwardMap.hpp"
#include "LinearOperator.hpp"
#include "../utils/Quadrature.hpp"
#include "../utils/BrentsMethod.hpp"


template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
ForwardMap<Np, Formulation, Ps, Pang, Nroot>::ForwardMap(double delta, ClosedCurve& curve, 
    double wavelengths_per_patch, double patch_split_wavenumber, int near_singular_patch_est, double p) :
        delta_(delta), curve_(curve), near_singular_patch_est_(near_singular_patch_est), p_(p)
{   
    // Precompute Quadrature nodes and weights for non-singular integration
    fejer_nodes_ = Eigen::ArrayXd(Np); fejer_weights_ = Eigen::ArrayXd(Np);
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

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
void ForwardMap<Np, Formulation, Ps, Pang, Nroot>::init_points_and_patches(ClosedCurve& curve, 
    double wavelengths_per_patch, double wavenumber) {

    std::vector<double> patch_lims = curve.compute_patch_lims(wavelengths_per_patch, wavenumber);
    
    num_patches_ = patch_lims.size() - 1;
    total_num_unknowns_ = num_patches_ * Np;

    patches_.clear(); patches_.reserve(num_patches_);

    xs_.clear(); ys_.clear(); nxs_.clear(); nys_.clear();
    xs_.reserve(num_patches_ * total_num_unknowns_); 
    ys_.reserve(num_patches_ * total_num_unknowns_);
    nxs_.reserve(num_patches_ * total_num_unknowns_); 
    nys_.reserve(num_patches_ * total_num_unknowns_); 

    for (long long ii = 0; ii < num_patches_; ii++) {
        patches_.emplace_back(patch_lims[ii], patch_lims[ii + 1], curve, delta_);

        patches_[ii].comp_bounding_box();

        patches_[ii].point_t_vals_.reserve(Np);
        patches_[ii].curve_jac = Eigen::ArrayXd(Np);

        // Store t, x, y, and curve jacobian for each point in the patch
        for (int jj = 0; jj < Np; jj++) { 
            double t_curr = ab2cd(fejer_nodes_[jj], -1.0, 1.0, patch_lims[ii], patch_lims[ii + 1]);

            patches_[ii].point_t_vals_.push_back(t_curr);

            Eigen::Vector2d normal = curve_.normal_t(t_curr);

            double normal_norm = normal.norm();

            patches_[ii].curve_jac(jj) = normal_norm;

            normal = normal / normal_norm;

            xs_.push_back(curve.xt(t_curr));
            ys_.push_back(curve.yt(t_curr));
            nxs_.push_back(normal(0));
            nys_.push_back(normal(1));
        }
    }
}

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
void ForwardMap<Np, Formulation, Ps, Pang, Nroot>::determine_patch_near_singular_points() {
    for (long long ii = 0; ii < num_patches_; ii++) {
    // Pointer to the current patch to test the near singular points.
    Patch<Nroot>* p_tp = &patches_[ii];
        for (long long jj = ii - near_singular_patch_est_; jj <= ii + near_singular_patch_est_; jj++) {

            long long test_patch_ind;

            // If there are only 2 patches this prevents double counting the patch
            if (num_patches_ <= 2 && jj > ii) { continue;}

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
            long long tp_val_start = test_patch_ind * Np;

            for (long long kk = tp_val_start; kk < tp_val_start + Np; kk++) {

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
                        p_tp->near_singular_point_lookup_.insert(kk);
                        p_tp->near_singular_point_ts_.push_back(p_tp->t1);

                    } else if (f(p_tp->t2) < delta_) {

                        p_tp->near_singular_point_indices_.push_back(kk);
                        p_tp->near_singular_point_lookup_.insert(kk);
                        p_tp->near_singular_point_ts_.push_back(p_tp->t2);

                    } else {

                        std::array<double,2> out = brent_find_minima(f, p_tp->t1, p_tp->t2);
             
                        // Check if the point in the bounding box is near singular.
                        if (out[1] < delta_) {
                            p_tp->near_singular_point_indices_.push_back(kk);
                            p_tp->near_singular_point_lookup_.insert(kk);
                            p_tp->near_singular_point_ts_.push_back(out[0]);
                        }

                    }
                    
                }

            }

        }

    }
}

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
Eigen::VectorXcd ForwardMap<Np, Formulation, Ps, Pang, Nroot>::single_patch_point_compute_precomputations
    (const Patch<Nroot>& patch, double tsing, double xsing, double ysing,
        std::complex<double> wave_number, Eigen::Ref<Eigen::VectorXcd> out_precomps) {

    // Set up vector of Cheb left and Cheb right
    Eigen::ArrayXd cheb_nm1  =  Eigen::ArrayXd::Constant(num_ns, 1.0);  
    Eigen::ArrayXd cheb_n(num_ns);    Eigen::ArrayXd cheb_right_n(num_ns);
    Eigen::ArrayXd cheb_np1(num_ns);  Eigen::ArrayXd cheb_right_np1(num_ns);

    // Store the change of variables, and its derivative evaluated at the cheb points
    Eigen::ArrayXd xi_alpha(num_ns), dxi_alpha(num_ns);
    
    // Change of variablesfrom [t1,t2] to [-1,1] cov and curve jacobians. 
    Eigen::ArrayXd g_cov(num_ns), curve_jac(num_ns);
    
    // Set up Green function evaluations on the left and right
    Eigen::ArrayXcd  green(num_ns);
    
    // Location in [-1,1] of the singularity
    double alpha = ab2cd(tsing, patch.t1, patch.t2, -1.0, 1.0);

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
            xi_alpha(ii)  = xi_cov_center(fejer_nodes_ns_(ii), alpha, p_);
            dxi_alpha(ii) = dxi_cov_center(fejer_nodes_ns_(ii), alpha, p_);
        } else if (sing_loc == -1) {
            xi_alpha(ii)  = xi_cov_left(fejer_nodes_ns_(ii), p_);
            dxi_alpha(ii) = dxi_cov_left(fejer_nodes_ns_(ii), p_);
        } else {
            xi_alpha(ii)  = xi_cov_right(fejer_nodes_ns_(ii), p_);
            dxi_alpha(ii) = dxi_cov_right(fejer_nodes_ns_(ii), p_);
        }

        cheb_n(ii) = xi_alpha(ii);

        // Compute the change of variables from [-1,1] to t1, t2;
        double g_xi_alpha = ab2cd(xi_alpha(ii), -1.0, 1.0, patch.t1, patch.t2);

        const double xt = curve_.xt(g_xi_alpha);
        const double yt = curve_.yt(g_xi_alpha);
  
        // Compute the curve Jacobian
        Eigen::Vector2d normal = curve_.normal_t(g_xi_alpha);
        curve_jac(ii) = normal.norm();

        // Compute the Green's Function
        normal       /= curve_jac(ii);
        green(ii) = GF<Formulation>(xsing, ysing, xt, yt, normal(0), normal(1),
        eta_, wave_number);
    }

    Eigen::ArrayXcd integrand_wo_cheb = 
        ((patch.t2 - patch.t1) / (2.0)) * green * curve_jac * dxi_alpha * fejer_weights_ns_;

    // No need to element wise multiply by all ones.
    out_precomps(0) = integrand_wo_cheb.sum();

    out_precomps(1) = (integrand_wo_cheb * cheb_n).sum();

    for (int ii = 2; ii < Np; ii++) {
        cheb_np1 = (2.0 * xi_alpha  * cheb_n) - cheb_nm1;

        out_precomps(ii) = (integrand_wo_cheb * cheb_np1).sum();

        // Update the Chebyshev polynomials
        cheb_nm1.swap(cheb_n);   
        cheb_n.swap(cheb_np1);
    }

    return out_precomps;

}

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
void ForwardMap<Np, Formulation, Ps, Pang, Nroot>::compute_precomputations(std::complex<double> wavenumber) {
    // Loop over patches
    for (long long ii = 0; ii < num_patches_; ii++) {
        Patch<Nroot>& patch = patches_[ii];

        patch.precomputations_ = 
            Eigen::MatrixXcd(Np, Np + patch.near_singular_point_ts_.size());
        
        // Singular Precomputations
        for (int jj = 0; jj < Np; jj++) {
            double xsing = xs_[ii * Np + jj];
            double ysing = ys_[ii * Np + jj];

            single_patch_point_compute_precomputations(patch, patch.point_t_vals_[jj],
            xsing, ysing, wavenumber, patch.precomputations_.col(jj));
        }

        // Near Singular Precomputations
        for (size_t jj = 0; jj < patch.near_singular_point_indices_.size(); jj++) {

            double tsing = patch.near_singular_point_ts_[jj];
            double xsing = xs_[patch.near_singular_point_indices_[jj]];
            double ysing = ys_[patch.near_singular_point_indices_[jj]];
   

            single_patch_point_compute_precomputations(patch, tsing, xsing, ysing,
            wavenumber, patch.precomputations_.col(Np + jj));
        }

    }
    
}

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
void ForwardMap<Np, Formulation, Ps, Pang, Nroot>::compute_intensities(Eigen::VectorXcd& density) {
 
    for (long long ii = 0; ii < num_patches_; ii++) {
 
        double m121_jac = ab2cdjac(-1,1,patches_[ii].t1,patches_[ii].t2);

        density.segment<Np>(ii * Np).array() *= m121_jac * fejer_weights_ * patches_[ii].curve_jac;
    }

}


template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
Eigen::VectorXcd ForwardMap<Np, Formulation, Ps, Pang, Nroot>::
    compute_non_singular_interactions_unacc(Eigen::VectorXcd& density, std::complex<double> wave_number) {

    Eigen::VectorXcd out = Eigen::VectorXcd::Zero(Np * num_patches_);
    
    // Loop over all points by looping over patches (corr. to rows of A)
    // The outer loop goes over the rows, and the inner one the row vector times density product
    for (long long pind1 = 0; pind1 < num_patches_; pind1++) {
        long long xind_start = Np * pind1;
        for (long long xind = xind_start; xind < xind_start+Np; xind++) {
            double xtarg = xs_[xind]; 
            double ytarg = ys_[xind];

            // Loop over all patches (correspond to columns of A)
            for (long long pind2 = 0; pind2 < num_patches_; pind2++) {
                // Make sure it is not singular or near singular
                if ((pind1 != pind2) && (patches_[pind2].near_singular_point_lookup_.count(xind) == 0)) {
            
                    // Compute the Green function
                    Eigen::ArrayXcd green(Np);
    
                    for (int ii = 0; ii < Np; ii++) {
                        double xsource = xs_[Np * pind2 + ii];
                        double ysource = ys_[Np * pind2 + ii];
                        double xnsource = nxs_[Np * pind2 + ii];
                        double ynsource = nys_[Np * pind2 + ii];

                        green(ii) = GF<Formulation>(xtarg, ytarg, xsource, ysource, 
                            xnsource, ynsource, eta_, wave_number);

                    }
                    
                    out(xind) += (green * density.segment<Np>(Np * pind2).array()).sum();
                }
            }
        }
    }

    return out;
}


template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
Eigen::VectorXcd ForwardMap<Np, Formulation, Ps, Pang, Nroot>::compute_sing_near_sing_interactions
    (Eigen::VectorXcd& density) { 
    
    Eigen::VectorXcd out = Eigen::VectorXcd::Zero(Np * num_patches_);
    for (long long pind = 0; pind < num_patches_; pind++) {
        // Compute the Chebshev coefficients of the density corresponding to that patch

        Eigen::VectorXcd coeffs = Cheb1D::interp_1d<Np>(density.segment<Np>(pind * Np));
        
        // Multiply the Chebshev coefficients by the precomputations

        Eigen::VectorXcd precomps_sum =  coeffs.transpose() * patches_[pind].precomputations_;

        // Add to the singular parts
        out.segment<Np>(pind * Np).array() += precomps_sum.segment<Np>(0).array();

        // Add to the near singular parts
        for (size_t sing_ind = 0; sing_ind < patches_[pind].near_singular_point_indices_.size(); sing_ind++) {

            out(patches_[pind].near_singular_point_indices_[sing_ind]) += 
                precomps_sum(Np + sing_ind);
        }
    }

    return out;
}

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
Eigen::VectorXcd ForwardMap<Np, Formulation, Ps, Pang, Nroot>::compute_Ax_unacc
    (Eigen::VectorXcd& density, std::complex<double> wave_number) {
        
        // This must go before computing the intensities, as compute intensities includes the integration weights.
        Eigen::VectorXcd sns = compute_sing_near_sing_interactions(density);

        compute_intensities(density);

        Eigen::VectorXcd ns  = compute_non_singular_interactions_unacc(density, wave_number);

        sns += ns;
        return sns;
}

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
void ForwardMap<Np, Formulation, Ps, Pang, Nroot>::compute_Ax_unacc
(Eigen::VectorXcd& density, Eigen::VectorXcd& out, std::complex<double> wave_number) {
    out = compute_Ax_unacc(density, wave_number);
}

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
void ForwardMap<Np, Formulation, Ps, Pang, Nroot>::precomps_and_setup(std::complex<double> wavenumber, int nlevels, BoxTree<Ps,Pang>& boxes) {
    
    auto start = std::chrono::high_resolution_clock::now();
    compute_precomputations(wavenumber);
    auto end   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to compute precomputations: " << elapsed.count() << std::endl;

    boxes.Reinitialize(xs_, ys_, nxs_, nys_, nlevels, wavenumber);

    start = std::chrono::high_resolution_clock::now();
    boxes.CheckIfSingandNearSingInNeighborhood(xs_, ys_, patches_);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time taken to check sing and near sing in neighborhood: " << elapsed.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    init_sort_sing_point(boxes.getInverse());
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time taken to compute sort sing point: " << elapsed.count() << std::endl;

}

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
void ForwardMap<Np, Formulation, Ps, Pang, Nroot>::init_sort_sing_point( 
const std::vector<long long>& inverse) {
    sort_sing_point_.resize(total_num_unknowns_);


    for (int patch_idx = 0; patch_idx < num_patches_; patch_idx++) {
        // Build a single unordered_set for this patch
        std::unordered_set<long long> patch_set;

        int patch_point_start = patch_idx * Np;

        // Insert all points in the patch
        for (int i = patch_point_start; i < patch_point_start + Np; i++) {
            patch_set.insert(inverse[i]);
        }

        // Copy the same set to all points in this patch
        for (int i = patch_point_start; i < patch_point_start + Np; i++) {
            sort_sing_point_[inverse[i]] = patch_set;  // copy assignment
        }
    }

    // Near singular points. We track the points in that patch which each point is near singular to
    // which are computed in the precomputations
    for (int patch_idx = 0; patch_idx < num_patches_; patch_idx++) {

        int patch_point_start = patch_idx * Np;

        for (int i : patches_[patch_idx].near_singular_point_indices_) {
             
             // Add all of the points in the patch which the target point x is near singular to
             for (int j = patch_point_start; j < patch_point_start + Np; j++) {

                sort_sing_point_[inverse[i]].insert(inverse[j]);

             }
        }

    }
}

template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
Eigen::VectorXcd ForwardMap<Np, Formulation, Ps, Pang, Nroot>::compute_Ax_acc
    (Eigen::VectorXcd& density, BoxTree<Ps,Pang>& boxes) {
    
    auto start = std::chrono::high_resolution_clock::now();
    Eigen::VectorXcd sns = compute_sing_near_sing_interactions(density);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "time for compute sing near sing interactions " << elapsed.count() << std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    compute_intensities(density);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "time for compute intensities " << elapsed.count() << std::endl;

    // TODO: Speedup This makes a deep copy since vector always owns its data.
    // Could speed up by rewriting boxtree to take Eigen::VectorXcd
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::complex<double>> v_density(density.data(), density.data() + density.size());
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "time for deep copy of Eigen::vector to std;:vector " << elapsed.count() << std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    boxes.template Solve<Formulation>(true, v_density, sort_sing_point_);
    // boxes_.template Solve<Formulation>(true, v_density, sort_sing_point_);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "time for boxtree solve " << elapsed.count() << std::endl;

    Eigen::Map<Eigen::VectorXcd> ns(v_density.data(), v_density.size());

    sns += ns;

    return sns;
   
 }

 //////////////////////////////// GMRES Solution Functions///////////////////////////////
 template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
 Eigen::VectorXcd ForwardMap<Np, Formulation, Ps, Pang, Nroot>::solve_unacc
 (const Eigen::VectorXcd& rhs, std::complex<double> wavenumber) {
    
    long long N = rhs.size();

    compute_precomputations(wavenumber);

    auto A = [this, wavenumber](Eigen::VectorXcd& x, Eigen::VectorXcd& solution) -> 
    void {compute_Ax_unacc(x, solution, wavenumber);};

    LinearOperator B(N,N,A);

    Eigen::opGMRES<LinearOperator> func(B);

    func.setMaxIterations(GMRES_MAX_ITER_);
    func.set_restart(GMRES_MAX_ITER_);
    func.setTolerance(GMRES_TOLERANCE_);

    Eigen::VectorXcd x = func.solve(rhs);

    return x;
    
 }
//////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////// Propagation functions
 template <int Np, FormulationType Formulation, int Ps, int Pang, int Nroot>
 Eigen::VectorXcd ForwardMap<Np, Formulation, Ps, Pang, Nroot>::propagate_unacc
 (const std::vector<double> x_prop, const std::vector<double> y_prop, Eigen::VectorXcd density,
    std::complex<double> wavenumber) {

    compute_intensities(density);

    Eigen::VectorXcd out(x_prop.size());

    Eigen::ArrayXcd green_func_evals(xs_.size());

    for (size_t i = 0; i < x_prop.size(); i++) {

        out[i] = 0.0;

        for (size_t j = 0; j < xs_.size(); j++) {

            out[i] += density[j] * GF<Formulation>(x_prop[i], y_prop[i], xs_[j], ys_[j],
            nxs_[j], nys_[j], eta_, wavenumber);

        }

    }

    return out;
 }
