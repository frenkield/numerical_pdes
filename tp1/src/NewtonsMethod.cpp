#include <cmath>
#include <iostream>
#include "NewtonsMethod.hpp"

double NewtonsMethod::findNextPoint(double point) {
    return point - function(point) / derivative(point);
}

double NewtonsMethod::findRoot(double startPoint) {

    double error = 1;
    int iteration = 0;
    double root = startPoint;

    while (error > tolerance && iteration++ < maxIterations) {
        root = findNextPoint(root);
        error = fabs(function(root));
    }

    return root;
};
