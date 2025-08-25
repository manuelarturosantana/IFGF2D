#pragma once

#include <vector>
#include <cmath>
#include <functional>

//.Send x in ab to cd
double inline ab2cd(const double x, const double a, const double b, const double c, const double d)
{

    return ((d - c) / (b - a)) * (x - a) + c;

}

double inline ab2cdjac(const double a, const double b, const double c, const double d) {
    return ((d - c) / (b - a));
}

// Fejer 1st quadrature rule
void inline fejerquadrature1(std::vector<double>& nodes, std::vector<double>& weights, int N)
{

    for (int i = 0; i < N; i++) {

        nodes[i] = std::cos(M_PI * (2.0 * i + 1.0) / (2.0 * N));
        weights[i] = 0.0;

        for (int j = 1; j <= std::floor(N * 0.5); j++) {

            weights[i] += std::cos(j * M_PI * (2.0 * i + 1.0) / N) / (4.0 * j*j - 1.0);

        }

        weights[i] = (2.0 / N) * (1.0 - 2.0 * weights[i]);

    }

}

///////////////////////////// GSS ////////////////////////////////////////////////////////
double golden_section_search(const std::function<double (double)>& f, double a, double b, 
        double tolerance) {
    const double invphi = (std::sqrt(5.0) - 1) / 2.0;

    while (std::abs(b - a) > tolerance) {
        double c = b - (b - a) * invphi;
        double d = a + (b - a) * invphi;
        if (f(c) < f(d)) {
            b = d;
        } else {
           a = c;
        }
    }

    return (a + b) / 2.0;
}

/////////////////////////// RP Changes of Variables //////////////////////////////////////

/// v function used in defining the w change of variables.
///@param t input paramter
///@param p Strength of cancelation for the change of variables
double inline v_func(double t, double p) {
    double val = (0.5 - (1.0 / p)) * std::pow((t / M_PI) - 1.0,3);
    val += (1.0/p) * ((t / M_PI) - 1) + 0.5;
    return val;
}

// Derivative of v_func
double inline vp_func(double t, double p) {
    double val = (3.0 / M_PI) * (0.5 - (1.0 / p)) * std::pow((t / M_PI) - 1.0, 2);
    val += 1.0 / (p * M_PI);
    return val;
}

/// base form of w change of variables. t \in [0,2\pi]
double inline w_cov_base(double t, double p) {
    double val = 2.0 * M_PI * std::pow(v_func(t, p),p);
    val /= (std::pow(v_func(t,p),p) + std::pow(v_func(2*M_PI - t, p),p));
    return val;
}

// Derivative of the base form of the change of variables
double inline wp_cov_base(double t, double p) {
    double val = 2 * M_PI * p * vp_func(t,p);
    val *= std::pow(v_func(t,p), p - 1) * std::pow(v_func(2 * M_PI - t, p),p - 1);
    val /= std::pow(
        std::pow(v_func(t, p),p) + std::pow(v_func(2 * M_PI, p),p), 2
     );
     return val;
}

/// normalized form of the w change of variables. t in [-1,1]
/// (1/pi)w(pi t + pi)
double inline w_cov(double t, double p) {
    return M_1_PI * w_cov_base(M_PI * t + M_PI , p);
}

/// Derivative of the normalized form of the w change of variables
double inline wp_cov(double t, double p) {
    return wp_cov_base(M_PI * t + M_PI, p);
}


