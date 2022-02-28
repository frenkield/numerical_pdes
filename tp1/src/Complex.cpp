#include "Complex.hpp"

Complex operator+(const Complex& c1, const Complex& c2) {
    Complex sum;
    sum.real = c1.real + c2.real;
    sum.imaginary = c1.imaginary + c2.imaginary;
    return sum;
}

Complex operator-(const Complex& c1, const Complex& c2) {
    return c1 + (-1 * c2);
}

Complex operator*(const Complex& c1, const Complex& c2) {
    Complex product;
    product.real = (c1.real * c2.real) - (c1.imaginary * c2.imaginary);
    product.imaginary = (c1.real * c2.imaginary) + (c1.imaginary * c2.real);
    return product;
}

Complex operator*(double a, const Complex& c) {
    return Complex(c.real * a, c.imaginary * a);
}

Complex operator*(const Complex& c, double a) {
    return operator*(a, c);
}

Complex operator/(const Complex& c1, const Complex& c2) {
    Complex c2Bar = c2.conjugate();
    double divisor = 1 / (c2 * c2Bar).getReal();
    return (c1 * c2Bar) * divisor;
}

ostream& operator<<(ostream& stream, Complex const& complex) {

    stream.precision(4);

    if (abs(complex.getImaginary()) < 0.0001) {
        stream << complex.getReal() << " + 0i";

    } else {
        stream << complex.getReal() << " + " << complex.getImaginary() << "i";
    }

    return stream;
}
