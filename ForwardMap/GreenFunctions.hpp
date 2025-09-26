#pragma once
#include <complex>
#include <iostream>
#include <stdlib.h>

#include "../complex_bessel-master/include/complex_bessel.h"

/*
Things left TODO:
- Double layer close interactions - will probably require higher derivatives. Seems better 
- than doing interpolation, which gives you a lot to worry about for general curves
- Possibly the distinction between interior and exterior. It seemed to matter in sebastians code
*/

enum class FormulationType {
    SingleLayer,
    DoubleLayer,
    Combined,
    ThreeD
};

/// @brief Evaluates the kernel Green's function or it's normal derivative for the different layer
///         potentials. The current form is so it can be passed into the IFGF. Uses -ieta convention
/// @tparam Formulation: Which formulation to use defined by the enums above
/// @param p1x The x and y coordinates of the points.
/// @param p1y 
/// @param p2x 
/// @param p2y 
/// @param nx The x coordinate of the unit normal with respect to p2
/// @param ny The y coordiante of the unit normal with respect ot p2
/// @param coupling_parameter Couple parameter, only used if Combined is chosen
/// @param wavenumber Wave number
/// @return The value of the Green's function.

template<FormulationType Formulation>
std::complex<double> inline GF(const double p1x, const double p1y,
        const double p2x, const double p2y,
        const double nx, const double ny,
        const double coupling_parameter, std::complex<double> wavenumber)
{

    std::complex<double> solution;

    const double norm_diff = std::sqrt((p2x - p1x)*(p2x - p1x) + (p2y - p1y)*(p2y - p1y));

    if (norm_diff < 1e-14) {

        // Solution is initialized to zero
        solution = std::complex<double>(0.0,0.0);

    } else {   

        // Static means the variables retain their values between function calls!!!
        static std::complex<double> tmpfactor = std::complex<double>(0.0,1.0)/4.0;
        static std::complex<double> c_unit (0.0,1.0);

        // constexpr means only the relevant branch is compiled
        if constexpr (Formulation == FormulationType::SingleLayer) {
            solution = tmpfactor * sp_bessel::hankelH1(0, wavenumber * norm_diff);

        } else if constexpr (Formulation == FormulationType::DoubleLayer) {

            // Compute beta = < r, n > / |r|^2
            
            const double rDotNorm = (p1x - p2x) * nx - (p1y - p2y) * ny;
            solution = tmpfactor * wavenumber * sp_bessel::hankelH1(1, wavenumber * norm_diff) 
                        * rDotNorm;
            // TODO: Decide how we are going to do close double layer interactions?
            // double beta = rDotNorm / (norm_diff * norm_diff);

            // if ((norm_diff < 1.0E-5) && (std::abs(rDotNorm) < 1.0E-6)) {

            //     beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, nx, ny, nz, dxdsx, dxdsy, dxdsz, dxdtx, dxdty, dxdtz, 
            //                                     dxdsdsx, dxdsdsy, dxdsdsz, dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

            // }            

        // } else if constexpr (Formulation == FormulationType::Combined){
        //     std::complex<double> solution_S = tmpfactor * sp_bessel::hankelH1(0, wavenumber * distance);

        //     // Compute beta = < r, n > / |r|^2
        // const double rDotNorm = (p1x - p2x) * nx - (p1y - p2y) * ny;
        //      std::complex<double> solution_D = tmpfactor * wavenumber * sp_bessel::hankelH1(1, wavenumber * norm_diff) 
        //                 * rDotNorm;
        //      solution = solution_D - c_unit * coupling_parameter * solution_S;

        } else if constexpr (Formulation == FormulationType::ThreeD) {
            solution = std::exp(c_unit * wavenumber * norm_diff) / (4.0 * M_PI * norm_diff);
        } else {
            std::cout << "GF unknown layer potential specification";
            std::exit(1);
        }

    }

    return solution;

}

// TODO SPEEDUP: Consider using cosine and sine to compute this
// As well as returning the inverse already computed 
std::complex<double> inline static factorization(const double distance, std::complex<double> wavenumber)
{
    // static constexpr std::complex<double> imag_unit (0,1.0);
    // return std::exp(imag_unit * wavenumber * distance) / distance;

    static constexpr std::complex<double> imag_unit (0,1.0);
    return std::exp(imag_unit * wavenumber * distance) / std::sqrt(distance);
} 
