#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Chebyshev.h"

enum InterpolatorScheme 
{
    Chebyshev = 0
};

template <InterpolatorScheme scheme>
class Interpolator
{

    private:

        static constexpr int Ps = 3;
        static constexpr int Pang = 5;

        const std::array<double, Ps> xs_;
        const std::array<double, Pang> xang_;

        const std::array<std::array<double, Ps>, Ps> Ts_; 
        const std::array<std::array<double, Pang>, Pang> Tang_;

        const int P_;

    public:

        constexpr Interpolator() : 
            xs_(Functions::SetupChebyshevPoints<Ps>()),
            xang_(Functions::SetupChebyshevPoints<Pang>()),
            Ts_(Functions::SetupMultipleTn<Ps>()),
            Tang_(Functions::SetupMultipleTn<Pang>()),
            P_(Ps*Pang)
        {
        }

        void GetInterpolationPoint(const int iter, double& x, double& y) const
        {
            x = xs_[iter % Ps];
            y = xang_[iter / Ps];
        }

        constexpr int GetNInterpPoints() const {return P_;}

        constexpr int GetInterpId(const int is, const int itheta) const
        {
            return itheta*Ps + is;
        }

        inline double Interpolate(const double x, const double y, double const * const vals_begin) const
        {
            std::array<double, Ps> TargetTs; 
            std::array<double, Pang> TargetTtheta;

            Functions::SetupMultipleTnArr<Ps>(x, TargetTs);
            Functions::SetupMultipleTnArr<Pang>(y, TargetTtheta);

            double result = 0.0;
            
            for (int sangiter1 = 0; sangiter1 < Pang; sangiter1++) {
                for (int sraditer = 0; sraditer < Ps; sraditer++) {
            
                    const long long liniter = sangiter1*Ps + sraditer;
                    const double cheb = TargetTtheta[sangiter1]*TargetTs[sraditer];                
            
                    result += cheb * vals_begin[liniter];
                }
            }

            return result;
        }

        void GenerateInterpolant(double* const arr) const
        {     

            std::array<std::array<double, Pang>, Ps> Cij;

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