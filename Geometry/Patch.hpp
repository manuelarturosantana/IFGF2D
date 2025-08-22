#pragma once

// #include <Eigen>
#include<array>

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

        BoundingBox(double Ax, double Ay, double Bx, double By,
                    double Cx, double Cy, double Dx, double Dy);

        bool is_inside(double xp, double yp);
};

// /// @brief 
// /// @tparam C A class which is derived from curve
// /// @tparam N The number of points used in the computation of the bounding box
// template <typename C, int N>
// class Patch {
//     static_assert(std::is_base_of<ClosedCurve, C>::value, "C must inherit from ClosedCurve");

//     public:
//         double t1, t2; // Limits in parameterization of the curve.
//         bool ist1open, ist2open;  // Tracks if patch has open endpts
//         bool ist1corner, ist2corner; // Tracks it patch has corner endpoints



//         BoundingBox bounding_box_; 
//         C curve_; // TODO Should this be a pointer to a curve?

//     Patch(double t1, double t2, C curve) : t1(t1), t2(t2), curve_(curve) {
//         bb_chebpts = Cheb1D::SetupChebyshevPoints<N>();
//     }

//     void comp_bounding_box() {
//         std::array<double, N> bb_chebpts = Cheb1D::SetupChebyshevPoints;




//         // evaluate x and y at the endpoints, as well as at any critical points
//         // compute the x min and max, and save as the bounding box points.
//     }

// };