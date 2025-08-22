#include "../../Geometry/Patch.hpp"
#include <iostream>

int main() {

    // Test with box centered 
    BoundingBox box(1.0, 0.0, 1.0, 1.0, 0.0,1.0, 0.0, 0.0);

    std::cout << "Is point (0.5, 0.5) inside the box? " 
              << (box.is_inside(0.5, 0.5) ? "Yes" : "No") << std::endl;
    std::cout << "Is point (1.5, 0.5) inside the box? " 
              << (box.is_inside(1.5, 0.5) ? "Yes" : "No") << std::endl;
    std::cout << "Is point (0.5, 1.5) inside the box? " 
              << (box.is_inside(0.5, 1.5) ? "Yes" : "No") << std::endl;
    std::cout << "Is point (0.99,0.001) inside the box? " 
              << (box.is_inside(0.99, 0.001) ? "Yes" : "No") << std::endl;
    std::cout << "Is point (1.01,0.001) inside the box? " 
              << (box.is_inside(1.01, 0.001) ? "Yes" : "No") << std::endl;
    std::cout << "Is point (1.0, 0.0) inside the box? " 
              << (box.is_inside(1.0, 0.0) ? "Yes" : "No") << std::endl;
    std::cout << "Is point (1.0, 1.0) inside the box? " 
              << (box.is_inside(1.0, 1.0) ? "Yes" : "No") << std::endl;

    return 0;
}