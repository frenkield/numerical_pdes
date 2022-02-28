#include <iostream>
#include <math.h>
#include "../src/Complex.hpp"
#include "../src/NewtonsMethodComplex.hpp"

using namespace std;

/*
    set style increment user
    set style line 1 lc rgb 'red'
    set style line 2 lc rgb 'blue'
    set style line 3 lc rgb 'green'

    set style data points
    plot 'points.txt' using 1:2:3 linecolor variable pt 5 ps 5 t ''
*/

bool equalsApproximately(double x, double y) {
    return fabs(x - y) <= 0.00001;
}

void generate() {

    Complex one(1, 0);

    auto function = [=](Complex z) -> Complex {return z * z * z - one;};
    auto derivative = [](Complex z) -> Complex {return 3 * z * z;};
    NewtonsMethodComplex nm(function, derivative);

    double minReal = -1;
    double maxReal = 1;
    double minImaginary = -1;
    double maxImaginary = 1;

    int pointCount = 100;

    double incrementReal = (maxReal - minReal) / pointCount;
    double incrementImaginary = (maxImaginary - minImaginary) / pointCount;

    for (int i = 0; i < pointCount + 1; i++) {

        for (int j = 0; j < pointCount + 1; j++) {

            Complex point(minReal + j * incrementReal, minImaginary + i * incrementImaginary);

            Complex root = nm.findRoot(point);
            int attractor = 1;

            if (equalsApproximately(root.getReal(), 1)) {
                attractor = 1;
            } else if (root.getImaginary() > 0) {
                attractor = 2;
            } else if (root.getImaginary() < 0) {
                attractor = 3;
            } else {
                cout << "error: " << root << endl;
            }

            cout << point.getReal() << " " << point.getImaginary() << " " << attractor << endl;

        }
    }
}

int main() {
    generate();
}
