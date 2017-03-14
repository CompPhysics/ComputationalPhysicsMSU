#ifndef Complex_H
#define Complex_H
//   various include statements and definitions
#include <iostream>          // Standard ANSI-C++ include files
#include <new>

//  My own Complex class
class Complex
{
private:
  double re, im; // real and imaginary part
public:
  Complex ();                              // Complex c;
  Complex (double re = 0.0, double im = 0.0); // Definition of a complex variable;
  Complex& operator= (const Complex& c); // c = a;   //  equate two complex variables
  ~Complex () {}                        // destructor
  double   Re () const;        // double real_part = a.Re();
  double   Im () const;        // double imag_part = a.Im();
  double   abs () const;       // double m = a.abs(); // modulus
  friend Complex operator+ (const Complex& a, const Complex& b);
  friend Complex  operator- (const Complex& a, const Complex& b);
  friend Complex operator* (const Complex& a, const Complex& b);
  friend Complex operator/ (const Complex& a, const Complex& b);
};

#endif
