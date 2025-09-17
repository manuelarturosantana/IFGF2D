// File declaring functions which create different right hand sides

#include<vector>
#include<complex>
#include<cmath>
#include<iomanip>

#include<Eigen/Dense>

#include "../IFGF_Source/PolarCoordinates.h"

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