
#include "Curve.hpp"
#include <iostream>


double Curve::comp_segement_length(double a, double b) {
    std::vector<double> xps, yps, ab_nodes(num_integration_points);

    for (int ii = 0; ii < num_integration_points; ii++) {
        ab_nodes[ii] = ab2cd(integration_nodes[ii], -1.0, 1.0, a, b);
    }

    gpt_vec(ab_nodes, xps, yps);

    double seg_len = 0.0;
    for (int ii = 0; ii < num_integration_points; ii++) {
        seg_len += std::sqrt(xps[ii] * xps[ii] + yps[ii] * yps[ii]) * integration_weights[ii];
    }
 
    return seg_len * ab2cdjac(-1.0, 1.0, a, b);
 }


void Curve::ref_seg_len(double a, double b, double size_lim, std::vector<double>& seg_endpts) {
    double len = comp_segement_length(a, b);

    if (len <= size_lim) {
        seg_endpts.push_back(b);
        return;
    }

    double mid = (a + b) / 2.0;
    ref_seg_len(a, mid, size_lim, seg_endpts);
    ref_seg_len(mid, b, size_lim, seg_endpts);

}

std::vector<double> Curve::compute_patch_lims(double wavelengths_per_patch,
    double wavenumber) {
    
    this->set_integration_points();
    double wavelength = 2 * M_PI / wavenumber;
    // Compute the length of the curve via integration
    double curve_len = comp_segement_length(tlim1, tlim2);

    // Split according to wavelengths_per_patch
    int num_wave_lens = std::ceil(curve_len / wavelength);
    int num_intervals = std::ceil(num_wave_lens / wavelengths_per_patch);

    std::vector<double> patch_lims; patch_lims.reserve(num_intervals);
    // Because ref_seg_len only returns right endpoint we initialize the first left endpoint
    patch_lims.push_back(tlim1);

    double t_len = tlim2 - tlim1;

    for (int ii = 0; ii < num_intervals; ii++) {
        double left_endpoint = tlim1 + t_len * ii / static_cast<double> (num_intervals);
        double right_endpoint = tlim1 + t_len * (ii + 1) / static_cast<double> (num_intervals);

        // DEBUG: Print the left and right endpoints
        // std::cout << "Left endpoint: " << left_endpoint << ", Right endpoint: " << right_endpoint << std::endl;

        ref_seg_len(left_endpoint, right_endpoint, wavelength * wavelengths_per_patch, patch_lims);
    } 

    this->integration_data_cleanup();

    return patch_lims;
}


void Curve::xt_vec(const std::vector<double>& ts, std::vector<double>& xs) const {
    xs.clear(); xs.reserve(ts.size());
    for (double tval : ts) {
        xs.push_back(xt(tval));
    }
}


void Curve::yt_vec(const std::vector<double>& ts, std::vector<double>& ys) const {
    ys.clear(); ys.reserve(ts.size());
    for (double tval : ts) {
        ys.push_back(yt(tval));
    }
}


void Curve::xpt_vec(const std::vector<double>& ts, std::vector<double>& xps) const {
    xps.clear(); xps.reserve(ts.size());
    for (double tval : ts) {
        xps.push_back(xpt(tval));
    }
}


void Curve::ypt_vec(const std::vector<double>& ts, std::vector<double>& yps) const {
    yps.clear(); yps.reserve(ts.size());
    for (double tval : ts) {
        yps.push_back(ypt(tval));
    }
}

/////////////////////////////////////// Derived Geometry Classes ///////////////////////////////////////
// Circle x and y coordinate functions
double Circle::xt(double t) const {
    return radius * std::cos(t) + centerx;
}

double Circle::yt(double t) const {
    return radius * std::sin(t) + centery;
}

// Circle x' and y' coordinate functions
double Circle::xpt(double t) const {
    return -radius * std::sin(t);
}

double Circle::ypt(double t) const {
    return radius * std::cos(t);
}

// Kite x and y coordinate functions
double Kite::xt(double t) const {
    return A * std::cos(t) + B * std::cos(2.0 * t) + C;
}

double Kite::yt(double t) const {
    return D * std::sin(t);
}

// Kite x' and y' coordinate functions
double Kite::xpt(double t) const {
    return -A * std::sin(t) - 2.0 * B * std::sin(2.0 * t);
}

double Kite::ypt(double t) const {
    return D * std::cos(t);
}

/////////////////////////////////////// Line /////////////////////////////////////////////

double Line::xt(double t) const {
    return (x2 * (t + 1.0) - x1 * (t - 1)) / 2.0;
}

double Line::yt(double t) const {
    return (y2 * (t + 1.0) - y1 * (t - 1)) / 2.0;
}

double Line::xpt([[maybe_unused]] double t) const {
    return (x2 - x1) / 2.0;
}

double Line::ypt([[maybe_unused]] double t) const {
    return (y2 - y1) / 2.0;
}