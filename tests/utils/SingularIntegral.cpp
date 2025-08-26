#include "../../utils/Quadrature.hpp"
#include <cmath>
#include <iostream>


// int main() {
//     double x_sing = 0.1;
//     auto f = [x_sing](double x) ->double {return std::log(std::abs(x - x_sing));};

//     double a = -2, b = 2;
//     int N = 20, p = 6;

//     double num_int = near_sing_int(f, a, b, x_sing, N, p);

//     double true_int = (b - x_sing) * std::log(std::abs(x_sing - b)) + (x_sing - a) * 
//         std::log(std::abs(x_sing - a)) - b + a;
    
//     std::cout << "The error is " << std::abs(num_int - true_int) << std::endl;
//     std::cout << "True int" << true_int << std::endl;

//     return 0;
// }