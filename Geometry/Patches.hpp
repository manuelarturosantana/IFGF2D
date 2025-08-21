#pragma once
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