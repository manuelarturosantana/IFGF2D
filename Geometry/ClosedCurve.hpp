#pragma once
#include <vector>
#include <cmath>

#include "../utils/quadrature.hpp"

class ClosedCurve {
    public:
        // Returns x and y coordinates of points in space.
        virtual void xt(const std::vector<double>& ts, std::vector<double>& xs, std::vector<double>& ys) const = 0;
        // Returns tangential derivatives of points in space.
        virtual void xpt(const std::vector<double>& ts, std::vector<double>& xps, std::vector<double>& yps) const = 0;

        // Computes the length of the integration segment corresponding to the parameterization
        // between a and b
        double comp_segement_length(double a, double b);

        // Returns a vector of right endpoints dividing ab in half (in parameter space) until
        //  the segment length is less than size_lim
        void ref_seg_len(double a, double b, double size_lim, std::vector<double>& seg_endpts);

        /// @brief For the 2pi parameterized curve compute the t lims to split it into patches
        /// @param wavelengths_per_patch Number of wavelengths before each patch. usually 1
        /// @param wavenumber The wavenumber for the problem
        /// @param num_int_points Number of points to use in the integration. May need to be larger for larger curves
        /// @return A vector where v[i],v[i+1] gives the t limits for patch i, i = 0,...,v.size() - 1;
        std::vector<double> compute_patch_lims(double wavelengths_per_patch, double wavenumber);

        virtual ~ClosedCurve() = default;

        int num_integration_points;
        std::vector<double> integration_nodes;
        std::vector<double> integration_weights;
   

        ClosedCurve(int num_integration_points = 200) : num_integration_points(num_integration_points) {
            // initialize the integration nodes and weights
            integration_nodes.resize(num_integration_points);
            integration_weights.resize(num_integration_points);
            fejerquadrature1(integration_nodes, integration_weights, num_integration_points);
        }
  
};

class Circle : public ClosedCurve { 
    public:
        double centerx;
        double centery;
        double radius;

        Circle(double centerx=0.0, double centery=0.0, double radius=1.0) :  ClosedCurve(),
         centerx(centerx), centery(centery), radius(radius) { }

        void xt(const std::vector<double>& ts, std::vector<double>& xs, std::vector<double>& ys) const override;
        void xpt(const std::vector<double>& ts, std::vector<double>& xps, std::vector<double>& yps) const override;
};


class Kite : public ClosedCurve { 
    public:
        // Kite is defined as (x(t), y(t)) with
        // x(t) = Acos(t) + Bcos(2t) + C
        // y(t) = Dsin(t)
        double A, B, C, D;

        // Default parameters correspond to the Colton and Kress kite
        Kite(double A=1.0, double B=0.65, double C=-0.65, double D=1.5) : ClosedCurve(),
         A(A), B(B), C(C), D(D) {};

        void xt(const std::vector<double>& ts, std::vector<double>& xs, std::vector<double>& ys) const override;
        void xpt(const std::vector<double>& ts, std::vector<double>& xps, std::vector<double>& yps) const override;
};


