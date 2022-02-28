#include <iostream>
#include <vector>
#include <cmath>
#include "../GC.hpp"
#include "../MatriceDuLaplacien.hpp"

using namespace std;

// plot "test1.txt" w l, sin(pi * x)

void printVector(vector<double>& vector) {
    for (auto i = vector.begin(); i != vector.end(); ++i) {
        cout << *i << " ";
    }
    cout << endl;
}

void printForGnuplot(vector<double>& vector1, vector<double>& vector2) {

    assert(vector1.size() == vector2.size());

    for (int i = 0; i < vector1.size(); i++) {
        cout << vector1[i] << " " << vector2[i] << endl;
    }
}

void generateSinRightSide(vector<double>& f, double step) {

    double piSquared = M_PI * M_PI;
    double stepSquared = step * step;

    for (int i = 0; i < f.size(); i++) {
        f[i] = stepSquared * piSquared * sin(M_PI * i * step);
    }
}

void generateXVector(vector<double>& xVector, double step) {
    for (int i = 0; i < xVector.size(); i++) {
        xVector[i] = i * step;
    }
}

int main(int argc, const char **argv) {

    int n = 50;
    double xMin = 0;
    double xMax = 1;
    double step = (xMax - xMin) / (n - 1);

    // avec conditions dirichlet homogene
    MatriceDuLaplacien A(n);

    vector<double> solution(n);
    vector<double> f(n);

    generateSinRightSide(f, step);

    // applique conditions dirichlet homogene
    f[0] = f[f.size() - 1] = 0;

    GradientConjugue(A, &f[0], &solution[0], n, 1e-8, 10);

    vector<double> x(n);
    generateXVector(x, step);

    printForGnuplot(x, solution);

    return 0;
}
