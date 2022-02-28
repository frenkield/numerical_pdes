#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>
#include <math.h>

using namespace std;

class Complex {

private:

    double real = 0;
    double imaginary = 0;

public:

    Complex() {}
    Complex(double real, double imaginary) : real(real), imaginary(imaginary) {}

    double getReal() const {
        return real;
    }

    double getImaginary() const {
        return imaginary;
    }

    Complex conjugate() const {
        return(Complex(real, -imaginary));
    }

    double norm() {
        Complex c = conjugate();
        Complex product = operator*(*this, c);
        return sqrt(product.real + product.imaginary);
    }

    friend Complex operator+(const Complex& c1, const Complex& c2);
    friend Complex operator-(const Complex& c1, const Complex& c2);
    friend Complex operator*(const Complex& c1, const Complex& c2);
    friend Complex operator*(double a, const Complex& c2);
    friend Complex operator*(const Complex& c2, double a);
    friend Complex operator/(const Complex& c1, const Complex& c2);

    friend ostream& operator<<(ostream& stream, Complex const& c);
};

#endif //COMPLEX_CPP
