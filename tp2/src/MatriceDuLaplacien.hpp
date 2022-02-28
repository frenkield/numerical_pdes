#ifndef MATRICEDULAPLACIEN_H
#define MATRICEDULAPLACIEN_H

#include <iostream>
#include <cassert>
#include <cmath>
#include "GC.hpp"

using namespace std;

class MatriceDuLaplacien : public MatVirt {

public:

    int n;
    bool isDirichletHomogeneous = true;

    MatriceDuLaplacien(int n) : MatVirt(n, n), n(n) {}

    friend std::ostream& operator<<(std::ostream& stream, MatriceDuLaplacien const& matrix);

    double operator()(int i, int j) const {

        assert(i < n && j < n);

        double value = 0;

        if (isDirichletHomogeneous && (i == 0 || i == n - 1)) {
            if (j == i) {
                value = 1;
            }

        } else if (i == j) {
            value = 2;

        } else if (fabs(i - j) == 1) {
            value = -1;
        }

        return value;
    }

    double *addmatmul(double *x, double *Ax) const {

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                Ax[i] += operator()(i, j) * x[j];
            }
        }

        return Ax;
    }
};

#endif