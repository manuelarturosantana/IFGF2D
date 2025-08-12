#include <cmath>
#include <iostream>

#include "BoxTree.h"

const double wavenumber_ = M_PI;

inline static void fct(const double x1, const double x2,
                       const double y1, const double y2,
                       const double densityreal, const double densityimag, 
                       double& phireal, double& phiimag) 
{

    static constexpr double tmpfactor = 1.0/M_PI/4.0;

    const double distance = std::sqrt((y1 - x1) * (y1 - x1) + (y2 - x2) * (y2 - x2));
    const double distance_inv = 1.0 / distance;

    const double kerreal = tmpfactor * cos(wavenumber_ * distance) * distance_inv;
    const double kerimag = tmpfactor * sin(wavenumber_ * distance) * distance_inv;

    phireal = densityreal * kerreal - densityimag * kerimag;
    phiimag = densityreal * kerimag + densityimag * kerreal;

}

inline static void fac(const double distance, double& re, double& im)
{

    re = cos(wavenumber_ * distance) / distance;
    im = sin(wavenumber_ * distance) / distance;

}

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

    for (long long iter = 0; iter < N; iter++) {
        
        L2error += (exact_real[iter] - approx_real[iter]) * (exact_real[iter] - approx_real[iter]) + 
                   (exact_imag[iter] - approx_imag[iter]) * (exact_imag[iter] - approx_imag[iter]);

        L2normsol += exact_real[iter]*exact_real[iter] + exact_imag[iter]*exact_imag[iter];

    }

    return std::sqrt(L2error/L2normsol);

}

int main(int argc, char *argv[]) 
{

    try {

        const long long N = 10000;
        const int nlevels = 6;

        const bool compute_singular_interactions = true;
        const bool compute_error = true;

        std::vector<double> pointsx(N);
        std::vector<double> pointsy(N);

        GenerateCircle(N, 1.0, pointsx, pointsy);    

        BoxTree boxes(pointsx, pointsy, nlevels, wavenumber_);

        std::vector<double> intensities_real(N, 1.0);
        std::vector<double> intensities_imag(N, 0.0);

        std::vector<double> tmpintensities_real(N, 1.0);
        std::vector<double> tmpintensities_imag(N, 0.0);
        
        boxes.Solve<&fct, &fac>(compute_singular_interactions, intensities_real, intensities_imag);

        if (compute_error) {    

            double error = ComputeError(intensities_real, intensities_imag, pointsx, pointsy, tmpintensities_real, tmpintensities_imag);

            std::cout << "L2 Error = " << error << std::endl;

        }
        
    } catch (std::exception& e) {

        std::cout << e.what() << std::endl;
        return 1;

    }

    return 0;

}