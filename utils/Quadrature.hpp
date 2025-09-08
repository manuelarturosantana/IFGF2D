#pragma once

#include <vector>
#include <cmath>
#include <functional>
#include <type_traits>


#include <iostream>
#include <Eigen/Dense>

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

template <typename T>
typename std::enable_if<std::is_same<T, Eigen::VectorXd>::value || std::is_same<T, Eigen::ArrayXd>::value, void>::type
 inline fejerquadrature1(T& nodes, T& weights, int N)
{

    for (int i = 0; i < N; i++) {

        nodes(i) = std::cos(M_PI * (2.0 * i + 1.0) / (2.0 * N));
        weights(i) = 0.0;

        for (int j = 1; j <= std::floor(N * 0.5); j++) {

            weights(i) += std::cos(j * M_PI * (2.0 * i + 1.0) / N) / (4.0 * j*j - 1.0);

        }

        weights(i) = (2.0 / N) * (1.0 - 2.0 * weights(i));

    }

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
    double vp = std::pow(v_func(t, p), p);
    double val = 2.0 * M_PI * vp;
    val /= (vp + std::pow(v_func(2 * M_PI - t, p),p));
    return val;
}

// Derivative of the base form of the change of variables
double inline wp_cov_base(double t, double p) {
    double val = 2.0 * M_PI * p * vp_func(t,p);
    val *= std::pow(v_func(t,p), p - 1.0) * std::pow(v_func(2 * M_PI - t, p),p - 1.0);
    val /= std::pow(
        std::pow(v_func(t, p), p) + std::pow(v_func(2 * M_PI - t, p),p), 2.0
     );
     return val;
}

// Taken from https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


/// @brief The xi change of variables used to cancel the singularities.
/// TODOSPEEDUP: If this doesn't get vectorized it may be too slow. Eigen could help vectorize this.
///              Also could combine the computation of the value and the change of variables to not recompute M_PI*std::abs(t) etc.
///              Althought I suspect that this latter speedup will be nearly unoticable.
/// @param t Value in [-1, 1]
/// @param alpha Value of the singularity in (-1,1)
/// @param p W change of variables parameter
/// @return The change of variables parameter.
double inline xi_cov_center(double t, double alpha, double p) {
    return alpha + ((sgn(t) - alpha) / M_PI) * w_cov_base(M_PI * std::abs(t), p);
} 

// alpha = 1
double inline xi_cov_right(double t, double p) {
    return 1.0 - ((2.0 / M_PI) * w_cov_base(M_PI * std::abs((t - 1.0) / 2.0),p));
}

// alpha = -1;
double inline xi_cov_left(double t,  double p) {
    return -1.0 + ((-2.0 / M_PI) * w_cov_base(M_PI * std::abs((t + 1.0) / 2.0),p));
}

//////////////////////////////// derivatives of the change of variables //////////////////
double inline dxi_cov_center(double t, double alpha, double p) {
    return (1.0 - alpha * sgn(t)) * wp_cov_base(M_PI * std::abs(t), p);
} 

// alpha = 1
double inline dxi_cov_right(double t, double p) {
    double temp = (t - 1.0) / 2.0;
    return -1.0 * wp_cov_base(M_PI * std::abs(temp), p) * sgn(temp);
}

// alpha = -1;
double inline dxi_cov_left(double t,  double p) {
    double temp = (t + 1.0) / 2.0;
    return -1.0 * wp_cov_base(M_PI * std::abs(temp),p) * sgn(temp);
}





/// normalized form of the w change of variables. t in [-1,1]
/// (1/pi)w(pi t + pi) Used in Sabhrants splitting code, but for now not used here.
// double inline w_cov(double t, double p) {
//     return M_1_PI * w_cov_base(M_PI * t + M_PI, p);
// }

// /// Derivative of the normalized form of the w change of variables
// double inline wp_cov(double t, double p) {
//     return wp_cov_base(M_PI * t + M_PI, p);
// }



/// @brief Prototype function for computing a near singular integral
/// @param f function to be integrated in the interval a b 
/// @param x_sing Location of the singularity in the interval
/// @param N The number of points to use on each side of the integration interval
////@param p change of variables parameter
/// @return The value of the integral

// double near_sing_int(const std::function<double (double)>& f, double a, double b, double x_sing, 
//         int N, double p) {

//     Eigen::VectorXd nodes(N), weights(N); 

//     fejerquadrature1(nodes, weights, N);

//     Eigen::VectorXd w_left(N), w_right(N), wp_left(N), wp_right(N), fs_left(N), fs_right(N);
//     for (int ii = 0; ii < N; ii++) {
//         const double x_left  = (nodes(ii) + 1.0) / 2.0;
//         const double x_right = (nodes(ii) - 1.0) / 2.0;

//         w_left(ii) =  w_cov(-x_left, p);  wp_left(ii) = wp_cov(-x_left, p);    
//         w_right(ii) = w_cov(x_right, p);  wp_right(ii) = wp_cov(x_right,p);

//         fs_left(ii) = f(x_sing - ((x_sing - a) * w_left(ii)));
//         fs_right(ii) = f(x_sing + ((b - x_sing) * w_right(ii)));
//     }

   
//     fs_left = ((x_sing - a) / 2.0) * fs_left;
//     fs_right = ((b - x_sing) / 2.0) * fs_right;

//     return (weights.array() * fs_left.array() * wp_left.array()).sum() 
//         + (weights.array() * fs_right.array() * wp_right.array()).sum();

// }



