#pragma once

// #include <Eigen>
#include <array>
#include <memory>
#include <Eigen/Dense>

#include "../utils/Chebyshev1d.hpp"
#include "../utils/Quadrature.hpp"
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

        // Matrix containing the points varying with column index, and tchebyshev polynomials
        // varying with row index. This convention makes it easier to grab a set of coefficients
        // as Eigen is column major
        Eigen::MatrixXcd precomputations_;

        BoundingBox bounding_box_;
        
        // Stores the t-value of points within the patch for precomputations.
        std::vector<double> point_t_vals_;

        // Stores the index over all curve points of the near singular point indices
        std::vector<long long> near_singular_point_indices_;
        // Location of tvals corresponding to the near singular points.
        std::vector<double>  near_singular_point_ts_; 

        // Storing the curve as a reference means we must guarantee that the curve passed 
        // in outlives all the patch objects that use it. If this is an issue
        // Consider using std::shared_ptr instead.
        ClosedCurve& curve_;

        // Hash table checking with keys at point index, and values as location of 
        // of start of precomputations in near singular precomputations
        // Vector or Eigen matrix of near singular precomputations

    Patch(double t1, double t2, ClosedCurve& curve, double delta) : 
        t1(t1), t2(t2), curve_(curve), delta(delta) {};

    /// @brief Computes a box that bounds the patch within a delta tolerance
    void comp_bounding_box();
};



#include "Patch.tpp"