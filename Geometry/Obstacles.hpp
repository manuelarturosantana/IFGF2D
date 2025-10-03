// This file contains declares functions to create various types of scattering obstacles
// Currently to do multiple objects just stitch together the current objects
#pragma once

#include <memory>
#include <vector>

#include "Curve.hpp"

/// @brief Construct a rectangle from line object
/// @param curves Output: A vector of curve objects
/// @param curve_touch_tracker Output: A vector of Junction structs recording how the individual curves touch
/// @param corner1x, corner1y, corner2x, corner2y Coordinates of the corners defining the rectangles
void make_rect(std::vector<std::unique_ptr<Curve>>& curves, 
    std::vector<std::vector<Junction>>& curve_touch_tracker,
    double corner1x = -0.5, double corner1y=0.5, double corner2x = 0.5, double corner2y = -0.5);