#include <iostream>
#include <math.h>
#include "../src/Complex.hpp"
#include "../src/NewtonsMethod.hpp"
#include "../src/NewtonsMethodComplex.hpp"

using namespace std;

bool equalsApproximately(double x, double y) {
    return fabs(x - y) <= 0.00001;
}

void test1() {
    Complex z;
    assert(z.getReal() == 0);
    assert(z.getImaginary() == 0);
}

void test2() {
    Complex z(2, 3);
    assert(z.getReal() == 2);
    assert(z.getImaginary() == 3);
}

void testAdd() {
    Complex z1(1, 2);
    Complex z2(4, 3.5);
    Complex sum = z1 + z2;

//    printf("%p\n", &sum);

    assert(sum.getReal() == 5);
    assert(sum.getImaginary() == 5.5);
}

void testProduct() {

    Complex z1(1, 2);
    Complex z2(4, 3.5);

    Complex product = z1 * z2;
    assert(product.getReal() == -3 && product.getImaginary() == 11.5);

    product = 4 * z1;
    assert(product.getReal() == 4 && product.getImaginary() == 8);
}

void testDivide() {

    Complex z1(1, 2);
    Complex z2(4, 3.5);

    Complex quotient = z1 / z2;

    assert(equalsApproximately(quotient.getReal(), 0.38938) &&
           equalsApproximately(quotient.getImaginary(), 0.159292));
}

void testNorm() {
    Complex z(-1, 2);
    assert(z.norm() == sqrt(5));
}

void testNewton1() {
    auto function = [](double x) -> double {return x * x;};
    auto derivative = [](double x) -> double {return 2 * x;};
    NewtonsMethod nm(function, derivative);
    assert(equalsApproximately(function(nm.findRoot(1)), 0));
}

void testNewtonComplex1() {

    Complex one(1, 0);

    auto function = [=](Complex z) -> Complex {return z * z * z - one;};
    auto derivative = [](Complex z) -> Complex {return 3 * z * z;};
    NewtonsMethodComplex nm(function, derivative);

    Complex root = nm.findRoot(Complex(200, 20));

    assert(equalsApproximately(root.getReal(), 1));
    assert(equalsApproximately(root.getImaginary(), 0));
}

int main() {
    test1();
    test2();
    testAdd();
    testProduct();
    testDivide();
    testNorm();
    testNewton1();
    testNewtonComplex1();
}
