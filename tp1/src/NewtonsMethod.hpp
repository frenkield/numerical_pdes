#ifndef TEST1_NEWTONSMETHOD_HPP
#define TEST1_NEWTONSMETHOD_HPP

#include <functional>

using namespace std;

class NewtonsMethod {

private:

    ::function<double(double)> function;
    ::function<double(double)> derivative;
    double tolerance = 0.00000001;
    int maxIterations = 100;

    double findNextPoint(double point);

public:

    NewtonsMethod(::function<double(double)> function, ::function<double(double)> derivative) {
        this->function = function;
        this->derivative = derivative;
    };

    double findRoot(double startPoint);
};


#endif //TEST1_NEWTONSMETHOD_HPP
