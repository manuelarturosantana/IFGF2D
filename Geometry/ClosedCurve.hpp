#pragma once

#include <vector>
#include <cmath>

#include "../utils/Quadrature.hpp"

// Waring: It is assumed that the parameterization is positively oriented
class ClosedCurve {
    public:
        // Returns x and y coordinates of points in space.
        void gt_vec(const std::vector<double>& ts, std::vector<double>& xs, std::vector<double>& ys) const {
            xt_vec(ts, xs); yt_vec(ts, ys);
         }
        // Returns tangential derivatives of points in space.
        void gpt_vec(const std::vector<double>& ts, std::vector<double>& xps, std::vector<double>& yps)  const{
            xpt_vec(ts, xps); ypt_vec(ts, yps);
        }

        // Return seperately the x and y coordinates of the points in space for vector
        void xt_vec(const std::vector<double>& ts, std::vector<double>& xs) const;
        void yt_vec(const std::vector<double>& ts, std::vector<double>& ys) const;
        void xpt_vec(const std::vector<double>& ts, std::vector<double>& xps) const;
        void ypt_vec(const std::vector<double>& ts, std::vector<double>& yps) const;

        // Single coordinate versions.
        // TODO SPEEDUP: Since these are virtual they go through a vtable lookup, and probably 
        // won't be inlined. New c++ compilers do devirtualization, but if code is too slow
        // one could implement this superclass in a functional framework ??
        virtual double inline xt(double t) const = 0;
        virtual double inline yt(double t) const = 0;
        virtual double inline xpt(double t) const = 0;
        virtual double inline ypt(double t) const = 0;
        

        // Computes the length of the integration segment corresponding to the parameterization
        // between a and b
        double comp_segement_length(double a, double b);

        // Returns a vector of right endpoints dividing ab in half (in parameter space) until
        //  the segment length is less than size_lim.
        void ref_seg_len(double a, double b, double size_lim, std::vector<double>& seg_endpts);

        /// @brief For the 2pi parameterized curve compute the t limits of patches of size less than
        ///        wavelengths_per_patch. Sets and cleans up the integration points 
        ///        used while doing this
        /// @param wavelengths_per_patch Number of wavelengths before each patch. usually 1
        /// @param wavenumber The wavenumber for the problem
        /// @return A vector where v[i],v[i+1] gives the t limits for patch i, i = 0,...,v.size() - 1;
        std::vector<double> compute_patch_lims(double wavelengths_per_patch, double wavenumber);

        virtual ~ClosedCurve() = default;

        int num_integration_points; // Number of integration points for computing patch size
        std::vector<double> integration_nodes;
        std::vector<double> integration_weights;

        void inline set_integration_points() {
            // initialize the integration nodes and weights
            integration_nodes.resize(num_integration_points);
            integration_weights.resize(num_integration_points);
            fejerquadrature1(integration_nodes, integration_weights, num_integration_points);
        }

        // Free the memory used by integration nodes and weights.
        // This will be useful so that each patch can have its own curve, with low memory.
        void integration_data_cleanup() {
            std::vector<double>().swap(integration_nodes);
            std::vector<double>().swap(integration_weights);
        }
   
        ClosedCurve(int num_integration_points = 200) : num_integration_points(num_integration_points) {
            set_integration_points();
        }  
};

class Circle : public ClosedCurve { 
    public:
        double centerx;
        double centery;
        double radius;

        Circle(double centerx=0.0, double centery=0.0, double radius=1.0) :  ClosedCurve(),
         centerx(centerx), centery(centery), radius(radius) { }

        double xt(double t) const override;
        double yt(double t) const override;
        double xpt(double t) const override;
        double ypt(double t) const override;

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

        double xt(double t) const override;
        double yt(double t) const override;
        double xpt(double t) const override;
        double ypt(double t) const override;
};


