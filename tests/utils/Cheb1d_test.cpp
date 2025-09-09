#include <iostream>
#include <math.h>

#include <Eigen/Dense>

#include "../../utils/Chebyshev1d.hpp"

// The chebysev approximation of the odd function (x exp(kx)) decays a little slower, and as expected
// as the even coefficients are zero, as all the even polynomials are even.
// https://math.stackexchange.com/questions/1736763/if-a-function-is-odd-even-then-its-best-polynomial-approximation-is-also-odd-ev

int main() {
  
    // Compute test set of points for exp(cos(kx))
    double k = 5;
    // Ugly, but usefule for the compile time code
    int N_test = 100;
    constexpr int num_points =  10;

      /////////////////// Eigen Version/////////////////////////////////////////////////
    // Eigen::VectorXd test_points(N_test);
    // Eigen::VectorXd fvals_test(N_test);
    // // auto f = [](double x) {return std::exp(x); };
    // auto f = [](double x) { return x * (x - (1.0/4.0) * (x - 1.0 / 2.0)); };
    // // auto f = [k] (double x) {return std::exp(std::cos(k * x)); };

    // // Fill test points with equally spaced values in the interval [-1,1]
    // for (int i = 0; i < N_test; i++) {
    //     test_points(i) = -1 + 2.0 * i / (N_test - 1);
    //     fvals_test(i) = f(test_points(i));
    // }

    // std::array<double, num_points> xs = Cheb1D::setup_chebyshev_points<num_points>();
    // Eigen::VectorXcd fvals(num_points);

    // for (int i = 0; i < num_points; i++) {
    //     fvals(i) = f(ab2cd(xs[i],-1,1,-1,1));
    // }

    // Eigen::VectorXcd coeffs = Cheb1D::interp_1d<num_points>(fvals);

    // // Print the coefficients
    // std::cout << "Coefficients for N = " << num_points << ": " << "\n";
    // for (int j = 0; j < num_points; j++) {
    //     std::cout << coeffs[j] << "\n";
    // }
    // std::cout << std::endl;


    // return 0;
/////////////////////////////////// Vector Version ///////////////////////////////////////
std::vector<double> test_points(N_test);
    std::vector<double> fvals_test(N_test);
    // auto f = [](double x) {return std::exp(x); };
    auto f = [](double x) { return x * (x - (1.0/4.0) * (x - 1.0 / 2.0)); };
    // auto f = [k] (double x) {return std::exp(std::cos(k * x)); };

    // Fill test points with equally spaced values in the interval [-1,1]
    for (int i = 0; i < N_test; i++) {
        test_points[i] = -1 + 2.0 * i / (N_test - 1);
        fvals_test[i] = f(test_points[i]);
    }

    std::vector<double> coeffs = Cheb1D::comp_f_cheb_coeffs<double, num_points>(f, -1, 1);

    // Print the coefficients
    std::cout << "Coefficients for N = " << num_points << ": " << "\n";
    for (int j = 0; j < num_points; j++) {
        std::cout << coeffs[j] << "\n";
    }
    std::cout << std::endl;

    
    // Evaluate the interpolation at the test points
    std::vector<double> interp_vals = Cheb1D::eval_1d_interp(coeffs, test_points);

    std::cout << "Interpolation values " << std::endl;
    for (int j = 0; j < N_test; j++) {
        std::cout << "Fvals test[i] " << fvals_test[j] << " interp_vals " << interp_vals[j] << std::endl;
    }



    // Compute the max error
    double max_error = 0.0;
    for (int j = 0; j < N_test; j++) {
        double error = std::abs(interp_vals[j] - fvals_test[j]);
        if (error > max_error) {
            max_error = error;
        }
    }
    
    // print the error
    std::cout << "N = " <<  num_points << ", Max Error = " << max_error << std::endl;

    return 0;

    
}



//////////////////////////////////////////////////////////////////////////////////////////
    // Print test points
    // for (int i = 0; i < N_test; i++) {
    //     std::cout << "Test point " << i << ": " << test_points[i] << std::endl;
    // }
   
    // loop over N values
    // std::array<double, num_points> xs = Cheb1D::setup_chebyshev_points<num_points>();

    // // Print the x points
    // for (int i = 0; i < num_points; i++) {
    //     std::cout << "xs[i] " << i << ": " << xs[i] << std::endl;
    // }

    
    // std::vector<double> fvals(num_points);

    // // Evaluate exp(cos(k x)) on the interval [-1,1] at cheb points
    // for (int j = 0; j < num_points; j++) {
    //     // fvals[j] = xs[j] * std::exp(std::cos(k * xs[j]));
    //     fvals[j] = f(xs[j]);
    //     // fvals[j] = exp(xs[j]);
    // }

    // // Print the function values
    // for (int j = 0; j < num_points; j++) {
    //     std::cout << "Function value at Chebyshev point " << j << ": " << fvals[j] << std::endl;
    // }

    // std::vector<double> coeffs = Cheb1D::interp_1d(fvals, xs);
