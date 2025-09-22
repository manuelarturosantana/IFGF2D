#include <cmath>
#include <iostream>
#include <complex>
#include <random>
#include <numeric>
#include <complex>
#include <chrono>


#include "BoxTree.h"

#include "../complex_bessel-master/include/complex_bessel.h"
#include "../utils/DebugUtils.hpp"




int num_f_evals = 0; // Check the number of function evaluations, which experiementally is about N log N :)
const double wavenumber_ = M_PI;
const std::complex<double> imag_unit (0,1.0);

// This is the Hankel function case
inline static void fct(const double x1, const double x2,
                       const double y1, const double y2,
                       const double densityreal, const double densityimag, 
                       double& phireal, double& phiimag) 
{

    // static constexpr double tmpfactor = 1.0/M_PI/4.0;

    const double distance = std::sqrt((y1 - x1) * (y1 - x1) + (y2 - x2) * (y2 - x2));
    std::complex<double> ker = sp_bessel::hankelH1(0, wavenumber_ * distance);
    const double kerreal = std::real(ker);
    const double kerimag = std::imag(ker);
    // const double distance_inv = 1.0 / distance;

    // const double kerreal = tmpfactor * cos(wavenumber_ * distance) * distance_inv;
    // const double kerimag = tmpfactor * sin(wavenumber_ * distance) * distance_inv;

    phireal = densityreal * kerreal - densityimag * kerimag;
    phiimag = densityreal * kerimag + densityimag * kerreal;
    num_f_evals += 1;
}

inline static void fac(const double distance, double& re, double& im)
{
    // Doing this gives a couple more digits of accuracy
    re = cos(wavenumber_ * distance) / std::sqrt(distance);
    im = sin(wavenumber_ * distance) / std::sqrt(distance);

    // re = cos(wavenumber_ * distance);
    // im = sin(wavenumber_ * distance);
    

} 


// inline static void fct(const double x1, const double x2,
//                        const double y1, const double y2,
//                        const double densityreal, const double densityimag, 
//                        double& phireal, double& phiimag) 
// {

//     static constexpr double tmpfactor = 1.0/M_PI/4.0;

//     const double distance = std::sqrt((y1 - x1) * (y1 - x1) + (y2 - x2) * (y2 - x2));
//     const double distance_inv = 1.0 / distance;

//     const double kerreal = tmpfactor * cos(wavenumber_ * distance) * distance_inv;
//     const double kerimag = tmpfactor * sin(wavenumber_ * distance) * distance_inv;

//     phireal = densityreal * kerreal - densityimag * kerimag;
//     phiimag = densityreal * kerimag + densityimag * kerreal;
//     num_f_evals += 1;
    
// }

// // In the paper this is known as the centered factor, but without the 1 /(4 PI)
// inline static void fac(const double distance, double& re, double& im)
// {

//     re = cos(wavenumber_ * distance) / distance;
//     im = sin(wavenumber_ * distance) / distance;
    

// } 

void GenerateCircle(const long long npoints, const double r, std::vector<double>& x, std::vector<double>& y) {

    for (long long i = 0; i < npoints; i++) {

        x[i] = r * cos(2 * M_PI / npoints * i);
        y[i] = r * sin(2 * M_PI / npoints * i);

    }

}

double ComputeError(const std::vector<double>& approx_real, 
                    const std::vector<double>& approx_imag, 
                    const std::vector<double>& x, const std::vector<double>& y,
                    const std::vector<double>& intensity_real, const std::vector<double>& intensity_imag)
{

    const long long N = x.size();

    std::vector<double> exact_real(N, 0.0);
    std::vector<double> exact_imag(N, 0.0);

    double L2error = 0.0;
    double L2normsol = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    for (long long target_iter = 0; target_iter < N; target_iter++) {

        double tmpreal, tmpimag;

        for (long long source_iter = 0; source_iter < N; source_iter++) {

            if (target_iter == source_iter)
                continue;
            
            fct(x[source_iter], y[source_iter],
                x[target_iter], y[target_iter],
                intensity_real[source_iter], intensity_imag[source_iter],
                tmpreal, tmpimag);

            exact_real[target_iter] += tmpreal;
            exact_imag[target_iter] += tmpimag;

        }

    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Direct computation time: " << elapsed.count() << " seconds\n";

    for (long long iter = 0; iter < N; iter++) {
        
        L2error += (exact_real[iter] - approx_real[iter]) * (exact_real[iter] - approx_real[iter]) + 
                   (exact_imag[iter] - approx_imag[iter]) * (exact_imag[iter] - approx_imag[iter]);

        L2normsol += exact_real[iter]*exact_real[iter] + exact_imag[iter]*exact_imag[iter];

    }

    return std::sqrt(L2error/L2normsol);

}

std::vector<double> random_vector(std::size_t N, double min = 0.0, double max = 1.0) {
    std::vector<double> vec(N);
    std::random_device rd;
    std::mt19937 gen(0);
    std::uniform_real_distribution<> dis(min, max);

    for (auto& x : vec) {
        x = dis(gen);
    }
    return vec;
}

int main() 
{

    try {

        const long long N = 1000;
        const int nlevels = 6;

        const bool compute_singular_interactions = true;
        const bool compute_error = true;

        std::vector<double> pointsx(N);
        std::vector<double> pointsy(N);

        GenerateCircle(N, 1.0, pointsx, pointsy);    

        auto setup_start = std::chrono::high_resolution_clock::now();
        BoxTree boxes(pointsx, pointsy, nlevels, wavenumber_);
        auto setup_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> setup_elapsed = setup_end - setup_start;
        std::cout << "Setup time: " << setup_elapsed.count() << " seconds\n";

        std::vector<double> intensities_real = random_vector(N, 0.0, 1.0);
        std::vector<double> intensities_imag(intensities_real);

        std::vector<double> tmpintensities_real(intensities_real);
        std::vector<double> tmpintensities_imag(intensities_imag);

        // std::vector<double> intensities_real(N, 1.0);
        // std::vector<double> intensities_imag(N, 1.0);

        // std::vector<double> tmpintensities_real(N, 1.0);
        // std::vector<double> tmpintensities_imag(N, 1.0);


        auto solve_start = std::chrono::high_resolution_clock::now();
        boxes.Solve<&fct, &fac>(compute_singular_interactions, intensities_real, intensities_imag);
        auto solve_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> solve_elapsed = solve_end - solve_start;
        std::cout << "Solve time: " << solve_elapsed.count() << " seconds\n";

        if (compute_error) {    

            double error = ComputeError(intensities_real, intensities_imag, pointsx, pointsy, tmpintensities_real, tmpintensities_imag);

            std::cout << "L2 Error = " << error << std::endl;

        }

        std::cout << "Number of F evals " << num_f_evals << std::endl;
        
    } catch (std::exception& e) {

        std::cout << e.what() << std::endl;
        return 1;

    }

    return 0;

}