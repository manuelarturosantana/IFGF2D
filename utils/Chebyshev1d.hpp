#pragma once

// Functions for 1D Tchebshev interpolation modified
// from some of the IFGF code, which is set up to do
// Interpolation in two variables.


#include <array>
#include <vector>
#include <cmath>
#include <utility>
#include <omp.h>
#include <functional>
#include <iostream>

#include <Eigen/Dense>

#include "quadrature.hpp"

namespace Cheb1D {

    template <int N>
    constexpr std::array<double, N> setup_chebyshev_points()
    {
        std::array<double, N> x{};

        #pragma unenroll        
        for (int i = 0; i < N; i++) {  
            
            x[N-1-i] = std::cos(M_PI*(2.0*i+1.0)/(2.0*N));
            
        }   

        return std::move(x); 

    }

    // Returns an array of SIZE corresponding to the Chebyshev evaluations at x
    template<std::size_t SIZE>
    inline void setup_multiple_Tn_arr(double x, std::array<double, SIZE>& ret) noexcept
    {            
        ret[0] = 1.0;
        ret[1] = x;

        #pragma unenroll
        for (int i = 2; i < SIZE; i++) {

            ret[i] = 2.0*x*ret[i-1] - ret[i-2];
        
        }

    }        

    // Returns an array of Chebyshev polynomials up to T_n at all N Tchebyshev points 
    // The first index is the polynomial degree, the second index is the point.
    template <int N>
    constexpr inline std::array<std::array<double, N>, N > setup_multiple_Tn() noexcept
    {
        const std::array<double, N> x = setup_chebyshev_points<N>();

        std::array<std::array<double, N>, N > polys{};

        #pragma unenroll        
        for (int i = 0; i < N; i++) {

            polys[0][i] = 1.0;
            polys[1][i] = x[i];
        
        }

        #pragma unenroll
        for (int n = 2; n < N; n++) {

            for (int i = 0; i < N; i++) {

                polys[n][i] = 2.0*x[i]*polys[n-1][i] - polys[n-2][i];
            
            }
        
        }

        return std::move(polys);
        
    }

    // TODO SPEEDUP: This computes the Chebyshev values for each time the intepolation is called.
    // We could use the setup multiple Tn to make this faster :)
    
    /// @brief Comptes the coefficients of the Chebyshev interpolation, where the coefficients
    /// function values are aready computed at the Chebyshev points. Assumes the open open grid :)
    /// @tparam T Any numeric type float, double, complex<float>, complex<double> etc.
    /// @tparam U Either a vector on an array holding the Tchebshev points
    /// @param fvals The function values evaluated at the chebxs
    /// @param chebxs The chebyshev points.
    /// @return A vector containg the coefficients
    template <typename T, typename U>
    std::vector<T> interp_1d(const std::vector<T>& fvals, const U& chebxs) {
        // Do the first two iterations so there is not branching in the loop.
        int N = fvals.size();
        std::vector<double> chebnm1(N), chebn(N);
        std::vector<T> coeffs(N);

        // Initialize the first two terms of the Tchebyshev recurrence
        for (int k = 0; k < N; k++) {
            chebnm1[k] = 1.0;
            chebn[k]   = chebxs[k];
        }

        // Loop over the x values to compute the first two coefficients
        for (int k = 0; k < N; k ++) {
            coeffs[0] += fvals[k];
            coeffs[1] += fvals[k] * chebxs[k];
        }
        coeffs[0] *= (1.0 / N);
        coeffs[1] *= (2.0 / N);

        // Compute the rest of the coefficients. 
        // i ranges across the coefficients, and b_i is always 2 here.
        for (int i = 2; i < N; i++) {
            coeffs[i] = 0.0;
            // k ranges across the x values
            for (int k = 0; k < N; k++) { 
                double val = 2.0 * chebxs[k] * chebn[k] - chebnm1[k];
                coeffs[i] +=  fvals[k] * (T) val;

                chebnm1[k] = chebn[k];
                chebn[k]   = val;
            }

            coeffs[i] *= (2.0 / N);
        }

        return coeffs;
    }

    /// @brief Evaluate the 1D tchebyshev interpolation at the points xs. This implementation
    /// assumes that the point xs are in the interval [-1,1].
    /// @tparam T float, double, complex<float>, complex<double>
    /// @tparam U Either a vector or an array
    /// @param coeffs Coefficients computed from Interp1d
    /// @param xs The points at which to evaluate the interpolation, assumes the xs are in [-1,1]
    /// @returns The output array evaluated at the interpolation points
    template<typename T, typename U>
    std::vector<T> eval_1d_interp(const std::vector<T>& coeffs, const U& xs) {
        std::vector<T> out(xs.size());

        for (size_t xind = 0; xind < xs.size(); xind++) {
            double chebnm1 = 1.0;
            double chebn   = xs[xind];
            // Add the first two values
            out[xind] += coeffs[0] * chebnm1 + coeffs[1] * chebn;



            // Evaluate the sum via the recursive formula
            for (size_t i = 2; i < coeffs.size(); i++) {
                double val = 2.0 * xs[xind] * chebn - chebnm1;

                out[xind] += coeffs[i] * val;

                chebnm1 = chebn;
                chebn  = val;
            }
        }

        return out;
    }

    /// @brief Evaluate the 1D tchebyshev interpolation at the points xs. This implementation
    /// assumes that the point xs are in the interval [a,b] and scales them to be in [-1,1].
    template<typename T, typename U>
    std::vector<T> eval_1d_interp(const std::vector<T>& coeffs, const U& xs, double a, double b) {
        std::vector<double> xs_scaled; xs_scaled.reserve(xs.size());
        for (size_t xind = 0; xind < xs.size(); xind++) {
            // Scale the points to be in [-1,1]
            xs_scaled[xind] = ab2cd(xs[xind],a,b,-1,1);
        }
        return eval_1d_interp(coeffs, xs_scaled);
    }

    /// @brief Computes the coefficients of a 1D function f from 0 to ... N - 1
    /// @tparam T Output type of f.
    /// @tparam N Number of points to use in the Chebyshev expansion
    /// @param f Scalar valued function
    /// @param a Left endpoint of interpolation interval
    /// @param b Right endpoint of interpolation interval
    /// @return A vector containing the coefficients of f
    template <typename T, int N>
    std::vector<T> comp_f_cheb_coeffs(const std::function<T (double)>& f, double a, double b) {
        std::array<double, N> xs = setup_chebyshev_points<N>();
        std::vector<T> fvals(N);

        for (int i = 0; i < N; i++) {
            fvals[i] = f(ab2cd(xs[i],-1,1,a,b));
        }

        return interp_1d(fvals, xs);

    }

    ///@brief Computes the roots of F via chebyshev interpolation
    template <int N>
    std::vector<double> cheb_root_finder(std::function<double (double)> f, double a, double b) {
        std::vector<double> f_coeffs = comp_f_cheb_coeffs<double, N>(f, a, b);

        // Compute the colleague matrix see chapter 18 of Trefethen approximation theory 
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N - 1, N - 1);
        A.diagonal(1).setConstant(0.5);   // superdiagonal
        A.diagonal(-1).setConstant(0.5);  // subdiagonal
        A(0,1) = 1.0;

        // TODO SPEEDUP: Eigen is column major, this is a last row traversal.
        // The matrices are small however, so many this won't be so bad.
        for (int i = 0; i < N - 1; i++) {
            //A is N-1 x N -1, eigen is zero indexed
            A(N-2, i) -= f_coeffs[i] / (2.0 * f_coeffs.back());
        }

        Eigen::EigenSolver<Eigen::MatrixXd> es(A);

        // Process Eigenvalues to only be those in the interval
        Eigen::VectorXcd eigenvalues = es.eigenvalues();
        std::vector<double> roots; roots.reserve(N - 1);

      
        for (int i = 0; i < eigenvalues.size(); i++) {
              // There will be complex roots, so we have to filter these out
            if(std::abs(eigenvalues(i).imag()) < 1e-3 &&
                eigenvalues(i).real() > -1 && eigenvalues(i).real() < 1) {

                roots.push_back(ab2cd(eigenvalues(i).real(),-1, 1, a, b));

            }
        }

        return roots;

    }
}