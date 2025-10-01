#include <iostream>
#include <functional>
#include <array>
#include <chrono>
#include<iomanip>

#include "../../complex_bessel-master/include/complex_bessel.h"

#include "../../ForwardMap/ForwardMap.hpp"
#include "../../Geometry/ClosedCurve.hpp"
#include "../../utils/Rhs.hpp"
#include "../../utils/Chebyshev1d.hpp"

using Clock = std::chrono::high_resolution_clock;


// Script to see if the runtime is scaling like nlogn
int main() {
    
    double delta = 0.01;
    double wave_lengths_per_patch = 1;
    int m = 3;
    constexpr int num_points = 10;


    std::vector<double> ks = {3 * 8, 3 * 16, 3 * 32, 3 * 64, 3 * 128, 3 * 256, 3 * 512, 3 * 1024};
    std::vector<double> npoints(ks.size());
    std::vector<double> total_times(ks.size());
    std::vector<double> IFGF_times(ks.size());

    for (size_t kind = 0; kind < ks.size(); kind++) {

        double k = ks[kind];

    // int nlevels = 4;

        Circle circle;
        
    
        // --- ForwardMap setup ---
        auto start_FM_setup = Clock::now();
        ForwardMap<num_points, FormulationType::SingleLayer, 5, 5> FM(delta, circle, wave_lengths_per_patch, k);
        auto end_FM_setup   = Clock::now();
        std::chrono::duration<double> dur_FM_setup = end_FM_setup - start_FM_setup;

        // Number of unknowns
        npoints[kind] = FM.total_num_unknowns_;
        // Cheecky rule of thumb with log base 4
        int nlevels = 3 + kind;
        // This is a hack to make it work for now
        if (kind == 6) nlevels -= 1;
        
        BoxTree<5,5> boxes;

        // --- Compute precomputations ---
        auto start_precomps = Clock::now();
        FM.compute_precomputations(k);
        auto end_precomps   = Clock::now();
        std::chrono::duration<double> dur_precomps = end_precomps - start_precomps;

        // --- Reinitialize boxes ---
        auto start_boxinit = Clock::now();
        boxes.Reinitialize(FM.xs_, FM.ys_, FM.nxs_, FM.nys_, nlevels, k);
        auto end_boxinit = Clock::now();
        std::chrono::duration<double> dur_box_init = end_boxinit - start_boxinit;

        // --- Check singular and near-singular ---
        auto start_check_sing = Clock::now();
        boxes.CheckIfSingandNearSingInNeighborhood(FM.xs_, FM.ys_, FM.patches_);
        auto end_check_sing   = Clock::now();
        std::chrono::duration<double> dur_check_sing = end_check_sing - start_check_sing;

        // --- Sort singular points ---
        auto start_sort_sing = Clock::now();
        FM.init_sort_sing_point(boxes.getInverse());
        auto end_sort_sing   = Clock::now();
        std::chrono::duration<double> dur_sort_sing = end_sort_sing - start_sort_sing;
        
        // --- RHS assembly ---
        Eigen::VectorXcd RHS = circle_eigenfunction(FM.xs_, FM.ys_, m);

        // --- Compute sing-near-sing interactions ---
        auto start_sns = Clock::now();
        Eigen::VectorXcd sns = FM.compute_sing_near_sing_interactions(RHS);
        auto end_sns   = Clock::now();
        std::chrono::duration<double> dur_sns = end_sns - start_sns;

        // --- Compute intensities ---
        auto start_intensities = Clock::now();
        FM.compute_intensities(RHS);
        auto end_intensities   = Clock::now();
        std::chrono::duration<double> dur_intensities = end_intensities - start_intensities;

        // --- Deep copy Eigen -> std::vector ---
        auto start_copy = Clock::now();
        std::vector<std::complex<double>> v_density(RHS.data(), RHS.data() + RHS.size());
        auto end_copy   = Clock::now();
        std::chrono::duration<double> dur_copy = end_copy - start_copy;

        // --- Boxtree solve ---
        auto start_solve = Clock::now();
        boxes.template Solve<FormulationType::SingleLayer>(true, v_density, FM.sort_sing_point_);
        auto end_solve   = Clock::now();
        std::chrono::duration<double> dur_solve = end_solve - start_solve;

        auto start_add = Clock::now();
        Eigen::Map<Eigen::VectorXcd> ns(v_density.data(), v_density.size());
        sns += ns;
        auto end_add   = Clock::now();
        std::chrono::duration<double> dur_add = end_add - start_add;

        // Print all durations
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "FM setup duration:           " << dur_FM_setup.count()      << " s\n";
        std::cout << "Precomputations duration:    " << dur_precomps.count()      << " s\n";
        std::cout << "Box init duration:           " << dur_box_init.count()      << " s\n";
        std::cout << "Check sing duration:         " << dur_check_sing.count()    << " s\n";
        std::cout << "Sort sing duration:          " << dur_sort_sing.count()     << " s\n";
        std::cout << "Sing-near-sing duration:     " << dur_sns.count()           << " s\n";
        std::cout << "Intensities duration:        " << dur_intensities.count()   << " s\n";
        std::cout << "Copy duration:               " << dur_copy.count()          << " s\n";
        std::cout << "BoxTree solve duration:      " << dur_solve.count()         << " s\n";
        std::cout << "Add duration:                " << dur_add.count()           << " s\n";
        std::cout << "---------------------------------------------" << std::endl;

        double IFGF_time = dur_box_init.count() +  dur_solve.count();
        double other_time = 
            dur_precomps.count() +
            dur_check_sing.count() +
            dur_sort_sing.count() +
            dur_sns.count() +
            dur_intensities.count() +
            dur_copy.count() +
            dur_add.count();
        
        total_times[kind] = other_time + IFGF_time;
        IFGF_times[kind]  = IFGF_time;
    
///////////////////////////////////////////////////////////////////////////////////////////////


//     std::complex<double> c_unit(0.0, 1.0);
//     std::complex<double> eval = ((c_unit * M_PI) / 2.0) * sp_bessel::besselJ(m,k) 
//                 * sp_bessel::hankelH1(m, k);
// ;
//     // The RHS gets changed :(
//     RHS = circle_eigenfunction(FM.xs_, FM.ys_, m);
   
//     Eigen::VectorXcd diff = Ax_acc - (eval * RHS);
//     std::cout << "Error in the Eigenfunction Test acc: " << diff.norm() << std::endl;

//     }
    }
    print_vector(total_times);
    print_vector(IFGF_times);
    print_vector(npoints);
    return 0;
}