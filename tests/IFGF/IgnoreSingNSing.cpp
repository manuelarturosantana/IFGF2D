#include <cmath>
#include <iostream>
#include <complex>
#include <random>
#include <numeric>
#include <complex>
#include <chrono>


#include "../../IFGF_Source/BoxTree.h"

#include "../../complex_bessel-master/include/complex_bessel.h"
#include "../../ForwardMap/GreenFunctions.hpp"
#include "../../utils/DebugUtils.hpp"

#include "../../ForwardMap/ForwardMap.hpp"
#include "../../Geometry/ClosedCurve.hpp"
#include "../../utils/Rhs.hpp"
#include "../../utils/Chebyshev1d.hpp"




int num_f_evals = 0; // Check the number of function evaluations, which experiementally is about N log N :)
const double wavenumber_ = M_PI;
const std::complex<double> imag_unit (0,1.0);

inline static std::complex<double> fct(const double x1, const double x2,
                       const double y1, const double y2,
                       const std::complex<double> density) 
{

    static constexpr double tmpfactor = 1.0/M_PI/4.0;
    static constexpr std::complex<double> imag_unit (0,1.0);

    const double distance = std::sqrt((y1 - x1) * (y1 - x1) + (y2 - x2) * (y2 - x2));
    const double distance_inv = 1.0 / distance;

    
    // const std::complex<double> ker = tmpfactor * std::exp(imag_unit * wavenumber_ * distance) * distance_inv;
    const double kerreal = tmpfactor * cos(wavenumber_ * distance) * distance_inv;
    const double kerimag = tmpfactor * sin(wavenumber_ * distance) * distance_inv;
    // This is slightly faster to do, but it won't matter in the Hankel function case.
    std::complex<double> ker = std::complex<double>(kerreal, kerimag);
    num_f_evals += 1;

       //     phireal = densityreal * kerreal - densityimag * kerimag;
//     phiimag = densityreal * kerimag + densityimag * kerreal;
    return density * ker;
 
    
}

// In the paper this is known as the centered factor, but without the 1 /(4 PI)
inline static std::complex<double> fac(const double distance)
{
    static constexpr std::complex<double> imag_unit (0,1.0);
    return std::exp(imag_unit * wavenumber_ * distance) / distance;

} 




template <FormulationType Formulation>
double ComputeError(const std::vector<std::complex<double>>& approx, 
    const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::complex<double>>& intensity, 
    std::vector<long long> inverse, std::vector<std::unordered_set<long long>> sort_sing_point)
{

    const long long N = x.size();

    std::vector<std::complex<double>> exact(N, 0.0);

    double L2error = 0.0;
    double L2normsol = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    for (long long target_iter = 0; target_iter < N; target_iter++) {

        // double tmpreal, tmpimag;

        for (long long source_iter = 0; source_iter < N; source_iter++) {

            if (target_iter == source_iter)
                continue;

            if (sort_sing_point[inverse[target_iter]].count(inverse[source_iter]) != 0) {
                // if (inverse[target_iter] == 118) {
                //     std::cout << "Comp_true_sol Skipping target " << inverse[target_iter] << " source " << inverse[source_iter] << std::endl;
                // }
                continue;
            }
            
            // For our text case normal derivative and coupling paramter are not used.
            exact[target_iter] += GF<Formulation>(x[source_iter], y[source_iter],
                x[target_iter], y[target_iter],
                0.0, 0.0, 1.0, wavenumber_) * intensity[source_iter];

        }

    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Direct computation time: " << elapsed.count() << " seconds\n";


    for (long long iter = 0; iter < N; iter++) {
        // Norm returns abs^2, which is what the code was doing
        L2error += std::norm(exact[iter] - approx[iter]);
        L2normsol += std::norm(exact[iter]);

    }

    return std::sqrt(L2error/L2normsol);

}
//////////////////////////// AI GENERATED UTILS//////////////////////////////////////////
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

std::vector<std::complex<double>> make_complex_vector(
    const std::vector<double>& re,
    const std::vector<double>& im)
{
    std::size_t n = re.size();
    std::vector<std::complex<double>> result;
    result.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        result.emplace_back(re[i], im[i]);
    }

    return result;
}
//////////////////////////////////////////////////////////////////////////////////////////

int main() 
{

    try {

        constexpr FormulationType ftype = FormulationType::SingleLayer;
        const int nlevels = 3;

        const bool compute_singular_interactions = true;
        const bool compute_error = true;

        double delta = 0.01;
        double k     = M_PI*200;
        // double k     = M_PI;
        double wave_lengths_per_patch = 1;

        Circle circle;
        constexpr int num_points = 10;

        ForwardMap<num_points, FormulationType::SingleLayer, 5, 5> FM(delta, circle, wave_lengths_per_patch, k);
        int N = FM.total_num_unknowns_;
        std::cout << "Total num unknowns " << N << std::endl;

        BoxTree<5,5> boxes_2;
        FM.precomps_and_setup(k, nlevels, boxes_2);


        std::vector<double> pointsx = FM.xs_;
        std::vector<double> pointsy = FM.ys_;  
        
        std::vector<double> pointsnx(pointsx);
        std::vector<double> pointsny(pointsy);

        BoxTree<10,10> boxes(pointsx, pointsy, pointsnx, pointsny, nlevels, wavenumber_);
      

        std::vector<double> intensities_real = random_vector(N, 0.0, 1.0);
        std::vector<double> intensities_imag(intensities_real);

        std::vector<std::complex<double>> intensities = make_complex_vector(intensities_real, intensities_imag);
        std::vector<std::complex<double>> tmpintensities(intensities);

      
        boxes.Solve<ftype>(compute_singular_interactions, intensities, FM.sort_sing_point_);
        // boxes.Solve<ftype>(compute_singular_interactions, intensities);


        if (compute_error) {    

            double error = ComputeError<ftype>(intensities, pointsx, pointsy, tmpintensities, boxes.getInverse(), FM.sort_sing_point_);

            std::cout << "L2 Error = " << error << std::endl;

        }

        std::cout << "Number of F evals " << num_f_evals << std::endl;

        // print_vector(boxes.getSorting());
        
    } catch (std::exception& e) {

        std::cout << e.what() << std::endl;
        return 1;

    }

    return 0;

}