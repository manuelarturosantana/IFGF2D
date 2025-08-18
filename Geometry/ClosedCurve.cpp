#include "ClosedCurve.hpp"
#include "../utils/quadrature.hpp"

 double ClosedCurve::comp_segement_length(double a, double b, int num_int_points) {
    std::vector<double> nodes(num_int_points), weights(num_int_points), xps, yps;
    fejerquadrature1(nodes, weights, num_int_points);

    for (int ii = 0; ii < num_int_points; ii++) {
        nodes[ii] = ab2cd(-1.0,1.0,a,b,nodes[ii]);
    }

    xpt(nodes, xps, yps);

    double seg_len = 0.0;
    for (int ii = 0; ii < num_int_points; ii++) {
        seg_len += std::sqrt(xps[ii] * xps[ii] + yps[ii] * yps[ii]) * weights[ii];
    }

    return seg_len;
 }

std::vector<double> ClosedCurve::compute_patch_lims(double wavelengths_per_patch,
    double wavenumber,int num_int_points = 100 ) {
    double wavelength = 2 * M_PI / wavenumber;
    // Compute the length of the curve via integration
   double curve_len = comp_segement_length(0, 2*M_PI, num_int_points);
    
   // Split according to wavelengths_per_patch
   int num_wave_lens = std::ceil(curve_len / wavelength);
   int num_intervals = std::ceil(num_wave_lens / wavelengths_per_patch);

   std::vector<double> patch_lims; patch_lims.reserve(num_intervals);
   
   // TODO: Finish this by putting in the patch intervals.
   for (int ii = 0; ii < num_intervals; ii++) {

   }
    
    // For each patch -> compute the length
    // If it is greater than wv_per_patch
    // divide the patch by two, then repeat
    // 
    
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
