#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <array>
#include <cmath>
#include <utility>

namespace Functions {

    template <int N>
    constexpr std::array<double, N> SetupChebyshevPoints()
    {
        std::array<double, N> x{};
              
        for (int i = 0; i < N; i++) {  
            
            x[N-1-i] = std::cos(M_PI*(2.0*i+1.0)/(2.0*N));
            
        }   

        return std::move(x); 

    }

    template<std::size_t SIZE>
    inline void SetupMultipleTnArr(double x, std::array<double, SIZE>& ret) noexcept
    {            
        ret[0] = 1.0;
        ret[1] = x;
        
        for (int i = 2; i < SIZE; i++) {

            ret[i] = 2.0*x*ret[i-1] - ret[i-2];
        
        }

    }        

    template <int N>
    constexpr inline std::array<std::array<double, N>, N > SetupMultipleTn() noexcept
    {
        const std::array<double, N> x = SetupChebyshevPoints<N>();

        std::array<std::array<double, N>, N > polys{};
               
        for (int i = 0; i < N; i++) {

            polys[0][i] = 1.0;
            polys[1][i] = x[i];
        
        }
        
        for (int n = 2; n < N; n++) {

            for (int i = 0; i < N; i++) {

                polys[n][i] = 2.0*x[i]*polys[n-1][i] - polys[n-2][i];
            
            }
        
        }

        return std::move(polys);
        
    }

}

#endif