#pragma once

// #include <Eigen>
#include<array>
#include <memory>

#include "../utils/Chebyshev1d.hpp"
#include "ClosedCurve.hpp"

// Class representing a box, with a method to check if a point is inside the box.
// 
//      __________
//  C  |          | B
//     |          |
//   D |_________ | A
// See https://stackoverflow.com/questions/2752725/finding-whether-a-point-lies-inside-a-rectangle-or-not
class BoundingBox {
    public:
        double Ax, Ay, Bx, By, Cx, Cy, Dx, Dy; // Coordinates of the corners
        double BmAx, BmAy, CmBx, CmBy; // Differences of vectors used for is_inside
        double BmAdotBmA, CmBdotCmB;   // squared norms of sides of the box

        BoundingBox() {};
        BoundingBox(double Ax, double Ay, double Bx, double By,
                    double Cx, double Cy, double Dx, double Dy);

        bool is_inside(double xp, double yp);
};

/// @brief 
/// @tparam N The number of points used root finding computation done while 
///         computing the bounding box.

template <int N>
class Patch {
    public:
        double t1, t2; // Limits in parameterization of the curve.
        double delta;  // RP near singularity parameter
        bool ist1open, ist2open;  // Tracks if patch has open endpts
        bool ist1corner, ist2corner; // Tracks it patch has corner endpoints

        BoundingBox bounding_box_; 
        std::unique_ptr<ClosedCurve> curve_; // TODO Should this be a pointer to a curve?

    Patch(double t1, double t2, std::unique_ptr<ClosedCurve> curve, double delta) : 
        t1(t1), t2(t2), curve_(std::move(curve)), delta(delta) {};

    void comp_bounding_box() {
        std::vector<double> xroots = Cheb1D::cheb_root_finder<N>([&](double t) { return curve_->xpt(t); }, t1, t2);
        std::vector<double> yroots = Cheb1D::cheb_root_finder<N>([&](double t) { return curve_->ypt(t); }, t1, t2);

        // TODO SPEEDUP: Avoid multiple calls to xt and yt
        double xmin = std::min(curve_->xt(t1), curve_->xt(t2));
        double xmax = std::max(curve_->xt(t1), curve_->xt(t2));
        double ymin = std::min(curve_->yt(t1), curve_->yt(t2));
        double ymax = std::max(curve_->yt(t1), curve_->yt(t2));

        for (double xr : xroots) {
            double xv = curve_->xt(xr);
            if (xv < xmin) xmin = xv;
            if (xv > xmax) xmax = xv;
        }

        for (double yr : yroots) {
            double yv = curve_->yt(yr);
            if (yv < ymin) ymin = yv;
            if (yv > ymax) ymax = yv;
        }

        xmin = xmin - delta; xmax = xmax + delta;
        ymin = ymin - delta; ymax = ymax + delta;

        // Ordering here follows the diagram for bounding box class
        bounding_box_ = BoundingBox(xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin);
    }

};