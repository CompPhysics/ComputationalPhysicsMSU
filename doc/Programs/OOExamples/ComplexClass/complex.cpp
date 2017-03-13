#include <iostream>
#include <cmath>
#include <complex>
using namespace std;
int main()
{
  std::complex<double> x(6.1,8.2), y(0.5,1.3);
  // write out x+y
  cout << x + y << x*y  << endl;
  return 0;
}
