#pragma once
#include "Patch.hpp"


// Helper function to subtract two vectors for the bounding box
void inline v1mv2(double x1, double y1, double x2, double y2, double& out1, double& out2) {
    out1 = x1 - x2;
    out2 = y1 - y2;
}

BoundingBox::BoundingBox(double Ax, double Ay, double Bx, double By,
            double Cx, double Cy, double Dx, double Dy) : Ax(Ax), Ay(Ay), Bx(Bx), By(By),
            Cx(Cx), Cy(Cy), Dx(Dx), Dy(Dy) {
    v1mv2(Bx, By, Ax, Ay, BmAx, BmAy);
    v1mv2(Cx, Cy, Bx, By, CmBx, CmBy);

    BmAdotBmA = BmAx * BmAx + BmAy * BmAy;
    CmBdotCmB = CmBx * CmBx + CmBy * CmBy;
}

bool BoundingBox::is_inside(double xp, double yp) {
    double pmAx, pmAy, pmBx, pmBy;
    v1mv2(xp, yp, Ax, Ay, pmAx, pmAy);
    v1mv2(xp, yp, Bx, By, pmBx, pmBy);

    double BmAdotpmA = BmAx * pmAx + BmAy * pmAy;
    double CmBdotpmB = CmBx * pmBx + CmBy * pmBy;

    return (0 <= BmAdotpmA && BmAdotpmA <= BmAdotBmA) &&
           (0 <= CmBdotpmB && CmBdotpmB <= CmBdotCmB);
}


/////////////////////////////////// Patch ////////////////////////////////////////////////

template <int N>
void Patch<N>::comp_bounding_box() {
    std::vector<double> xroots = Cheb1D::cheb_root_finder<N>([&](double t) { return curve_.xpt(t); }, t1, t2);
    std::vector<double> yroots = Cheb1D::cheb_root_finder<N>([&](double t) { return curve_.ypt(t); }, t1, t2);

    // TODO SPEEDUP: Avoid multiple calls to xt and yt
    double xmin = std::min(curve_.xt(t1), curve_.xt(t2));
    double xmax = std::max(curve_.xt(t1), curve_.xt(t2));
    double ymin = std::min(curve_.yt(t1), curve_.yt(t2));
    double ymax = std::max(curve_.yt(t1), curve_.yt(t2));

    for (double xr : xroots) {
        double xv = curve_.xt(xr);
        if (xv < xmin) xmin = xv;
        if (xv > xmax) xmax = xv;
    }

    for (double yr : yroots) {
        double yv = curve_.yt(yr);
        if (yv < ymin) ymin = yv;
        if (yv > ymax) ymax = yv;
    }

    xmin = xmin - delta; xmax = xmax + delta;
    ymin = ymin - delta; ymax = ymax + delta;

    // Ordering here follows the diagram for bounding box class
    bounding_box_ = BoundingBox(xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin);
}


