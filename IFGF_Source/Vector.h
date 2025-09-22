#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <algorithm>

/*
Assuming the vector A contains indices
And B contains values this function returns a sorted B,
and A returns incides containing the sort order.
*/

void new_sort(std::vector<long long>& vector_A, std::vector<long long>& vector_B) {

    stable_sort(vector_A.begin(), vector_A.end(),
        [&vector_B](long long i1, long long i2) {return vector_B[i1] < vector_B[i2];});

    std::vector<long long> tmp_vector_B(vector_B.size());

    for (size_t i = 0; i < vector_B.size(); i++) {

        tmp_vector_B[i] = vector_B[vector_A[i]];
                
    }

    vector_B = tmp_vector_B;

}

// return a sorted version of the vector by the indices in sorting
void arrange_in_sorted_order(std::vector<long long>& sorting, std::vector<double>& vector ) {
    
    std::vector<double> tmp_vector(vector.size());
    for (size_t i = 0; i < vector.size(); i++) {

        tmp_vector[i] = vector[sorting[i]];
                
    }
    vector = tmp_vector;

}

#endif