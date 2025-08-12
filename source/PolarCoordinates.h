#ifndef POLARCOORDINATES_H
#define POLARCOORDINATES_H

#include <cmath>

namespace Functions {

    inline double Atan2(const double y, const double x) noexcept {

        if (x > 0 && y >= 0)
            return atan(y/x);
        else if (x > 0 && y < 0)
            return atan(y/x) + 2.0*M_PI;
        else if (x < 0 && y > 0)
            return atan(y/x) + M_PI;
        else if (x < 0 && y == 0)
            return +M_PI;
        else if (x < 0 && y < 0)
            return atan(y/x) + M_PI;
        else if (x == 0 && y > 0)
            return M_PI_2;
        else if (x == 0 && y < 0)
            return 3.0*M_PI_2;
        else
            return 0;
            
    }

    inline void CartToPol(double& x, double& y) noexcept {

        const double r = std::sqrt(x*x + y*y);

        if (r == 0.0) {
            x = 0.0; y = 0.0;
            return;
        }

        const double theta = Atan2(y, x);

        x = r;
        y = theta;

    }

    inline void PolToCart(double & r, double & theta) noexcept {

        double sintheta, costheta;
        
        sincos(theta, &sintheta, &costheta);
        
        const double x = r * costheta;
        const double y = r * sintheta;

        r = x;
        theta = y;

    }

};

#endif