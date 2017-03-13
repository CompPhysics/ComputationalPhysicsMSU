#include <iostream>
#include <cmath>
#include "Complex.h"
using namespace std;
int main()
{
  Complex a(0.1,1.3);    // we declare a complex variable a
  Complex b(3.0), c(5.0,-2.3);  // we declare  complex variables b and c
  Complex d = b;         //  we declare  a new complex variable d
  d = a*c + b/a;  //   we add, multiply and divide two complex numbers
  cout << "Re(d)=" << d.Re() << ", Im(d)=" << d.Im() << endl;  // write out of the real and imaginary parts
  return 0;
}
