#include <memory>

#include "Obstacles.hpp"


void make_rect(std::vector<std::unique_ptr<Curve>>& curves, 
    std::vector<std::vector<Junction>>& curve_touch_tracker,
    double corner1x, double corner1y, double corner2x, double corner2y) 
    {   


        curves.emplace_back(std::make_unique<Line>(corner2x, corner2y, corner2x, corner1y));
        curves.emplace_back(std::make_unique<Line>(corner2x, corner1y, corner1x, corner1y));
        curves.emplace_back(std::make_unique<Line>(corner1x, corner1y, corner1x, corner2y));
        curves.emplace_back(std::make_unique<Line>(corner1x, corner2y, corner2x, corner2y));

        curve_touch_tracker = {
            {Junction(1, false), Junction(3, false)}, 
            {Junction(0, false), Junction(2, false)},
            {Junction(1, false), Junction(3, false)},
            {Junction(2, false), Junction(0, false)}
        };
        
    }