# IFGF README
This contains some Q&A from Sebastian on how the code works.

Q: What is the a non-relevant cone

A: There are 2 different names for the cones in the code (following Christoph's code notation)

Non-relevant cones are the cones obtained from a box setting. Initially, in the 2d case, you have 1x4 non-relevant cones (radius x theta, in polar coordinates)

Relevant cones are all the cones used from all the boxes at certain level. Initially, here you will have more than 4 (it depends on the number of boxes, cones that you use, etc.)

For example, at level D, for morton box 0 you can have 2 cones that you will use, for morton box 1 you can have 1 cone that you will use, etc. The relevant cone notation for them will be 0, 1, 2, etc; but the non-relevant cone notation will be between 0 and 3 (because you have 4 max)

The function from the picture (Car2Nonrelconefinds, given a box, the non-relevant cone where the point (x,y) is located

Q: What is protobox

A: That is used in the propagation
The idea is the following: given a box, you have certain cones that you will use (they are the keys of the 1st unordered map)

For each cone, you will have, for example, 3x5 interpolation points. For each interpolation point, you need to find the cones in the previous level that contain them

So, let's take interpolation point 0. First, let's get the children of the big box (that explains the std::array<, 4>)

For each child box, let's find the child cone that contains interpolation point 0

Suppose it is child cone number 3 (remember that these are non-relevant cones). Then, you go to the second unordered map, and for key 3 you add 0 to the vector

In summary, it is a data structure that contains the position of the interpolation points in the previous level

And you need to store this once since all the boxes will have the same scenario