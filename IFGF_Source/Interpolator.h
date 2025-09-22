#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include <complex>
#include "Chebyshev.h"

template <int Ps = 3, int Pang  =5>
class Interpolator
{

    private:

        const std::array<double, Ps> xs_;
        const std::array<double, Pang> xang_;

        const std::array<std::array<double, Ps>, Ps> Ts_; 
        const std::array<std::array<double, Pang>, Pang> Tang_;

        const int P_ = Ps * Pang; // The total number of interpolation points Ps * Pang;

    public:

        constexpr Interpolator() : 
            // This is some fancy templating stuff to make all of this evaluated at compile time.
            // TODO: For more flexible code we can call this precomputation cost. This is nice,
            // but means we have to come in an explicitly change everything, which is lame.
            xs_(Functions::SetupChebyshevPoints<Ps>()),
            xang_(Functions::SetupChebyshevPoints<Pang>()),
            Ts_(Functions::SetupMultipleTn<Ps>()),
            Tang_(Functions::SetupMultipleTn<Pang>()),
            P_(Ps*Pang)
        {
        }

        // 1D indexing to get the interpolation point. Uses leading dimension Ps, then each
        // row or column varyings the s variable.
        void GetInterpolationPoint(const int iter, double& x, double& y) const
        {
            x = xs_[iter % Ps];
            y = xang_[iter / Ps];
        }

        constexpr int GetNInterpPoints() const {return P_;}

        // Leading dimension is Ps
        constexpr int GetInterpId(const int is, const int itheta) const
        {
            return itheta*Ps + is;
        }

        //
        // double const * const means a constant pointer to a constant double
        // const at the end means that this function does not modify the objects member variables.
        // Evaluate the interpolant at some point x and y, (which I believe should really be s and theta).
        inline std::complex<double> Interpolate(const double x, const double y, std::complex<double> const * const vals_begin) const
        {
            std::array<double, Ps> TargetTs; 
            std::array<double, Pang> TargetTtheta;
            
            // Compute the Tchebyshev values for the x and y values.
            Functions::SetupMultipleTnArr<Ps>(x, TargetTs);
            Functions::SetupMultipleTnArr<Pang>(y, TargetTtheta);

            std::complex<double> result = 0.0;
            
            for (int sangiter1 = 0; sangiter1 < Pang; sangiter1++) {
                for (int sraditer = 0; sraditer < Ps; sraditer++) {
            
                    const long long liniter = sangiter1*Ps + sraditer;
                    const double cheb = TargetTtheta[sangiter1]*TargetTs[sraditer];                
            
                    result += cheb * vals_begin[liniter];
                }
            }

            return result;
        }

        // Compute the coefficients in the Chebyshev interpolant.
        // store them in array
        void GenerateInterpolant(std::complex<double>* const arr) const
        {     

            std::array<std::array<std::complex<double>, Pang>, Ps> Cij;

            for (int m = 0; m < Ps; m++) {
                for (int n = 0; n < Pang; n++) {

                    Cij[m][n] = 0.0;

                    for (int k = 0; k < Ps; k++) {
                        for (int l = 0; l < Pang; l++) {

                            Cij[m][n] += arr[GetInterpId(k, l)] * Ts_[m][k] * Tang_[n][l];

                        }
                    }

                }
            }

            for (int m = 0; m < Ps; m++) {
                for (int n = 0; n < Pang; n++) {

                    const int id = GetInterpId(m, n);

                    arr[id] = Cij[m][n];

                    double alphan = 2.0;
                    double alpham = 2.0;

                    if (n == 0) alphan = 1.0;
                    if (m == 0) alpham = 1.0;

                    arr[id] *= alpham * alphan / (Ps * Pang);

                }
            }

        }

};

#endif