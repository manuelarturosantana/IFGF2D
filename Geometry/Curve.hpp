#pragma once

#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include "../utils/Quadrature.hpp"

// Struct to track, along with a vector which curves touch which
struct Junction {
    int touching_curve; // Index of the other curve which this one touches.
    bool t1_touches_t1; // if true, the t1_lim of this patch touches the t1 lim of the other curve.
                        // otherwise the t1_lim of this patch the t2_lim of the other patch

    Junction(int touching_curve, bool t1_touches_t1) : touching_curve(touching_curve), 
        t1_touches_t1(t1_touches_t1) {};
};


// Waring: It is assumed that the parameterization is positively oriented
class Curve {
    public:

        // The user will have to take care to to change is_closed, tlim1, or tlim2 after a curve
        // has been initialized. However making them run time constants allows for the closed curves
        // to automatically be used as open curves
        bool is_closed;
        double tlim1, tlim2;

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

        virtual ~Curve() = default;

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
   
        Curve(bool is_closed, double tlim1, double tlim2, int num_integration_points = 200) : 
        is_closed(is_closed), tlim1(tlim1), tlim2(tlim2), num_integration_points(num_integration_points) {
            set_integration_points();
        }
        
        // Return an eigen vector of the normal (not unit), assuming a counter clockwise oriented curve.
        Eigen::Vector2d inline normal_t(double t) {
            Eigen::Vector2d normal;
            normal(0) = ypt(t);
            normal(1) = -xpt(t);
            return normal;
        };
        
};

class Circle : public Curve { 
    public:
        double centerx;
        double centery;
        double radius;

        Circle(bool is_closed=true, double tlim1 = 0, double tlim2 = 2.0 *M_PI, 
            double centerx=0.0, double centery=0.0, double radius=1.0) :  
            Curve(is_closed, tlim1, tlim2),
         centerx(centerx), centery(centery), radius(radius) { 
            tlim1 = 0;
            tlim2 = 2* M_PI;
         }

        double xt(double t) const override;
        double yt(double t) const override;
        double xpt(double t) const override;
        double ypt(double t) const override;

};


class Kite : public Curve { 
    public:
        // Kite is defined as (x(t), y(t)) with
        // x(t) = Acos(t) + Bcos(2t) + C
        // y(t) = Dsin(t)
        double A, B, C, D;

        // Default parameters correspond to the Colton and Kress kite
        Kite(bool is_closed=true, double tlim1 = 0, double tlim2 = 2.0 *M_PI, 
             double A=1.0, double B=0.65, double C=-0.65, double D=1.5) : 
             Curve(is_closed, tlim1, tlim2),
            A(A), B(B), C(C), D(D) {
         };

        double xt(double t) const override;
        double yt(double t) const override;
        double xpt(double t) const override;
        double ypt(double t) const override;
};

class Line : public Curve {
    public :
        // Line connection the points (x1,y1), (x2,y2)
        double x1, y1, x2, y2;

        Line(double x1, double y1, double x2, double y2) : Curve(false, -1.0, 1.0),
            x1(x1), y1(y1), x2(x2), y2(y2) {};

            double xt(double t) const override;
            double yt(double t) const override;
            double xpt(double t) const override;
            double ypt(double t) const override;

};

