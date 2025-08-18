#pragma once

#include <vector>
#include <cmath>

//.Send x in ab to cd
double inline ab2cd(const double a, const double b, const double c, const double d, const double x)
{

    return ((d - c) / (b - a)) * (x - a) + c;

}

// Fejer 1st quadrature rule
void fejerquadrature1(std::vector<double>& nodes, std::vector<double>& weights, int N)
{

    for (int i = 0; i < N; i++) {

        nodes[i] = std::cos(M_PI * (2.0 * i + 1.0) / (2.0 * N));
        weights[i] = 0.0;

        for (int j = 1; j <= std::floor(N * 0.5); j++) {

            weights[i] += std::cos(j * M_PI * (2.0 * i + 1.0) / N) / (4.0 * j*j - 1.0);

        }

        weights[i] = (2.0 / N) * (1.0 - 2.0 * weights[i]);

    }

}