#include <cmath>
#include "Complex.h"

Complex:: Complex () { re = im = 0.0; }

Complex:: Complex (double re_a, double im_a) {re = re_a; im = im_a; }

double Complex:: Re () const { return re; } //  getting the real part
double Complex:: Im () const { return im; }  //   and the imaginary part
double Complex:: abs () const { return sqrt(re*re + im*im); }
Complex& Complex:: operator= (const Complex& c)
{
  re = c.re;
  im = c.im;
  return *this;
}
inline Complex operator+ (const Complex& a, const Complex& b) { return Complex (a.re + b.re, a.im + b.im); }
inline Complex operator- (const Complex& a, const Complex& b) { return Complex (a.re - b.re, a.im - b.im); }

inline Complex operator* (const Complex& a, const Complex& b) {return Complex(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im);}
inline Complex operator/ (const Complex& a, const Complex& b) {return Complex( (a.re*b.re +a.im*b.im)/(b.re*b.re+b.im*b.im), (-a.im*b.re + a.im*b.im)/(b.re*b.re+b.im*b.im));}
