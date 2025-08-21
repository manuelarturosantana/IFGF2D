#include "Patches.hpp"


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