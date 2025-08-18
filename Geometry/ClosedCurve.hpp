#pragma once
#include <vector>
#include <cmath>

class ClosedCurve {
    public:
        // Returns x and y coordinates of points in space.
        virtual void xt(const std::vector<double>& ts, std::vector<double>& xs, std::vector<double>& ys) const = 0;
        // Returns tangential derivatives of points in space.
        virtual void xpt(const std::vector<double>& ts, std::vector<double>& xps, std::vector<double>& yps) const = 0;

        /// @brief For the 2pi parameterized curve compute the t lims to split it into patches
        /// @param wavelengths_per_patch Number of wavelengths before each patch. usually 1
        /// @param wavenumber The wavenumber for the problem
        /// @param num_int_points Number of points to use in the integration. May need to be larger for larger curves
        /// @return A vector where v[i],v[i+1] gives the t limits for patch i, i = 0,...,v.size() - 1;
        std::vector<double> compute_patch_lims(double wavelengths_per_patch, double wavenumber, int num_int_points = 100);

        // Computes the length of the integration segment corresponding to the parameterization
        // between a and b
        double comp_segement_length(double a, double b, int num_int_points);


        virtual ~ClosedCurve() = default;
};

class Circle : ClosedCurve { 
    public:
        double centerx;
        double centery;
        double radius;

        
        Circle(double centerx=0.0, double centery=0.0, double radius=1.0) : centerx(centerx),
            centery(centery), radius(radius) {}

        void xt(const std::vector<double>& ts, std::vector<double>& xs, std::vector<double>& ys) const override;
        void xpt(const std::vector<double>& ts, std::vector<double>& xps, std::vector<double>& yps) const override;

};


