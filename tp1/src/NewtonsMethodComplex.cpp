#include <cmath>
#include <iostream>
#include "NewtonsMethodComplex.hpp"

Complex NewtonsMethodComplex::findNextPoint(Complex point) {
    return point - function(point) / derivative(point);
}

Complex NewtonsMethodComplex::findRoot(Complex startPoint) {

    double error = 1;
    int iteration = 0;
    Complex root = startPoint;

    while (error > tolerance && iteration++ < maxIterations) {
        root = findNextPoint(root);
        error = function(root).norm();
    }

    return root;
};
