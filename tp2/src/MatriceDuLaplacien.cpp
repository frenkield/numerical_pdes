#include <iostream>
#include "MatriceDuLaplacien.hpp"

using namespace std;

std::ostream& operator<<(std::ostream& stream, MatriceDuLaplacien const& matrix) {

    for(int i = 0; i < matrix.n; i++) {
        for (int j = 0; j < matrix.m; j++) {
            stream << matrix(i, j) << " ";
        }
        stream << std::endl;
    }

    return stream;
}
