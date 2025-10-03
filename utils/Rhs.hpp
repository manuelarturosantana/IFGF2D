// File declaring functions which create different right hand sides

#include<vector>
#include<complex>
#include<cmath>
#include<iomanip>

#include<Eigen/Dense>

#include "../IFGF_Source/PolarCoordinates.h"
#include "../complex_bessel-master/include/complex_bessel.h"

/// @brief Creates the eigenfunction for the SLP evaluated on the circle. exp^(i m theta)
///        The eigenvalue is (i pi / 2)J_n(k)H_n^(1)(k);
/// @param xs Vector of x values to evaluate the eigenfunction on
/// @param ys vector of y values to evaluate the eigenfunction on
/// @param wave_number The wave number
/// @param m Parameter for the eigenfunction
/// @return A vector with the eigenvalues.
Eigen::VectorXcd circle_eigenfunction(const std::vector<double>& xs, 
    const std::vector<double> & ys, int m) {
    
        std::complex<double> c_unit(0.0,1.0);
        Eigen::VectorXcd out(xs.size());

        for (size_t ii = 0; ii < xs.size(); ii++) {
            double x = xs[ii]; double y = ys[ii];
            
            Functions::CartToPol(x, y);

            out(ii) = std::exp(c_unit * y * static_cast<double>(m));
        }

        return out;

}

/// @brief Compute the RHS given by a point source (2D Helmholtz green function)
/// @param xs x values on the boundary to evalute the RHS at
/// @param ys y values on the boundary to evalute the RHS at
/// @param x_0 x value of location of the point source
/// @param y_0 y value of the location of the point source
/// @param wavenumber Wave number of the fundamental solution
/// @return A vector containing the RHS evaluated at the wavenumber
Eigen::VectorXcd point_source(const std::vector<double>& xs, 
    const std::vector<double> & ys, double x_0, double y_0, std::complex<double> wavenumber) {

        Eigen::VectorXcd out(xs.size());

        for (size_t ii = 0; ii < xs.size(); ii++) {
            double x = xs[ii]; double y = ys[ii];
            
            double distance = std::sqrt(
                std::pow(x - x_0, 2) + std::pow(y - y_0,2)
            );

            out[ii] =  sp_bessel::hankelH1(0, wavenumber * distance);
        }

        return out;


    }

// Test case from Section 4.1 of The Helmholtz Dirichlet and Neumann problems on piecewise smooth open curves
// Helsing and Jiang
Eigen::VectorXcd test_case_straight_line(const std::vector<double>& xs, 
    const std::vector<double> & ys) 
    {

        Eigen::VectorXcd out(xs.size());

        for (size_t ii = 0; ii < xs.size(); ii++) {
            double x = xs[ii];

            out[ii] = 4.0 * std::pow(x,3)  + 2.0 * std::pow(x,2) - 3.0 * x - 1;
        }

        return out;


}