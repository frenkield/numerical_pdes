#ifndef NEWTONSMETHODCOMPLEX_HPP
#define NEWTONSMETHODCOMPLEX_HPP

#include <functional>
#include "Complex.hpp"

using namespace std;

class NewtonsMethodComplex {

private:

    ::function<Complex(Complex)> function;
    ::function<Complex(Complex)> derivative;
    double tolerance = 0.00000001;
    int maxIterations = 100;

    Complex findNextPoint(Complex point);

public:

    NewtonsMethodComplex(::function<Complex(Complex)> function, ::function<Complex(Complex)> derivative) {
        this->function = function;
        this->derivative = derivative;
    };

    Complex findRoot(Complex startPoint);
};


#endif
