#pragma once

#include <vector>

template <typename T>
void print_vector(const std::vector<T>& vec) {
    for (const auto& val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}