#include "ClosedCurve.hpp"
#include <iostream>
double ClosedCurve::comp_segement_length(double a, double b) {
    std::vector<double> xps, yps, ab_nodes(num_integration_points);

    for (int ii = 0; ii < num_integration_points; ii++) {
        ab_nodes[ii] = ab2cd(-1.0, 1.0, a, b, integration_nodes[ii]);
    }

    xpt(ab_nodes, xps, yps);

    double seg_len = 0.0;
    for (int ii = 0; ii < num_integration_points; ii++) {
        seg_len += std::sqrt(xps[ii] * xps[ii] + yps[ii] * yps[ii]) * integration_weights[ii];
    }
 
    return seg_len * ab2cdjac(-1.0, 1.0, a, b);
 }

void ClosedCurve::ref_seg_len(double a, double b, double size_lim, std::vector<double>& seg_endpts) {
    double len = comp_segement_length(a, b);

    if (len <= size_lim) {
        seg_endpts.push_back(b);
        return;
    }

    double mid = (a + b) / 2.0;
    ref_seg_len(a, mid, size_lim, seg_endpts);
    ref_seg_len(mid, b, size_lim, seg_endpts);

}


std::vector<double> ClosedCurve::compute_patch_lims(double wavelengths_per_patch,
    double wavenumber) {
    double wavelength = 2 * M_PI / wavenumber;
    // Compute the length of the curve via integration
    double curve_len = comp_segement_length(0, 2*M_PI);

    // Split according to wavelengths_per_patch
    int num_wave_lens = std::ceil(curve_len / wavelength);
    int num_intervals = std::ceil(num_wave_lens / wavelengths_per_patch);

    std::vector<double> patch_lims; patch_lims.reserve(num_intervals);
    // Because ref_seg_len only returns right endpoint we initialize the first left endpoint
    patch_lims.push_back(0.0);

    for (int ii = 0; ii < num_intervals; ii++) {
        double left_endpoint = 2 * M_PI * ii / static_cast<double> (num_intervals);
        double right_endpoint = 2 * M_PI * (ii + 1) / static_cast<double> (num_intervals);

        // DEBUG: Print the left and right endpoints
        std::cout << "Left endpoint: " << left_endpoint << ", Right endpoint: " << right_endpoint << std::endl;

        ref_seg_len(left_endpoint, right_endpoint, wavelength * wavelengths_per_patch, patch_lims);
    } 

    return patch_lims;
}



void Circle::xt(const std::vector<double>& ts, std::vector<double>& xs, std::vector<double>& ys) const {
    xs.clear(); ys.clear(); xs.reserve(ts.size()); ys.reserve(ts.size());
    for (double tval : ts) {
        xs.push_back(radius * std::cos(tval) + centerx);
        ys.push_back(radius * std::sin(tval) + centery);
    }
}

void Circle::xpt(const std::vector<double>& ts, std::vector<double>& xps, std::vector<double>& yps) const {
xps.clear(); yps.clear(); xps.reserve(ts.size()); yps.reserve(ts.size());
    for (double tval : ts) {
        xps.push_back(-radius * std::sin(tval) + centerx);
        yps.push_back(radius * std::cos(tval) + centery);
    }
}



void Kite::xt(const std::vector<double>& ts, std::vector<double>& xs, std::vector<double>& ys) const {
    xs.clear(); ys.clear(); xs.reserve(ts.size()); ys.reserve(ts.size());
    for (double tval : ts) {
        xs.push_back(A * std::cos(tval) + B * std::cos(2 * tval) + C);
        ys.push_back(D * std::sin(tval));
    }
}

void Kite::xpt(const std::vector<double>& ts, std::vector<double>& xps, std::vector<double>& yps) const {
    xps.clear(); yps.clear(); xps.reserve(ts.size()); yps.reserve(ts.size());
    for (double tval : ts) {
        xps.push_back(-A * std::sin(tval) - 2 * B * std::sin(2 * tval));
        yps.push_back(D * std::cos(tval));
    }
}
