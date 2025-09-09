#include <iostream>
#include <functional>
#include <array>
#include <chrono>


#include "../../ForwardMap/ForwardMap.hpp"
#include "../../Geometry/ClosedCurve.hpp"

// Script to Test the accuracy of the forward map on an eigenfunction test.

int main() {
    double delta = 0.1;
    double k     = 100;
    double wave_lengths_per_patch = 1.0;

    Kite kite;
    constexpr int num_points = 10;
    ForwardMap<num_points, FormulationType::SingleLayer, num_points> FM(delta, kite, wave_lengths_per_patch, k);

    auto start = std::chrono::high_resolution_clock::now();
    FM.compute_precomputations(std::complex<double>(k));
    auto end   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to compute precomputations: " << elapsed.count() << std::endl;


}