#pragma once
#include "ForwardMap.hpp"
#include "../utils/Quadrature.hpp"


template <int Np, int Nroot>
ForwardMap<Np, Nroot>::ForwardMap(double delta, ClosedCurve& curve, double wavelengths_per_patch, double patch_split_wavenumber) :
        delta(delta)
{

    initialze_points_and_patches(curve, wavelengths_per_patch, patch_split_wavenumber);
}

template <int Np, int Nroot>
std::vector<Patch<Nroot>> ForwardMap<Np, Nroot>::init_points_and_patches(ClosedCurve& curve, 
    double wavelengths_per_patch, double wavenumber) {
    
    std::vector<double> patch_lims = curve.compute_patch_lims(wavelengths_per_patch, wavenumber);

    num_patches_ = patch_lims.size() - 1;
    total_num_unknowns = num_patches_ * Np;

    patches_.clear() patches_.reserve(num_patches_);
    xs_.clear(); ys.clear(); 
    
    xs.reserve(num_patches_ * total_num_unknowns); 
    ys.reserve(num_patches_ * total_num_unknows);

    for (size_t ii = 0; ii < num_patches_; ii++) {
        patches_.emplace_back(patch_lims[ii], patch_lims[ii + 1], curve, delta_);
        patches_[i].comp_bounding_box();

        // Add the points to the patch
    }
}

