#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
// use namespace for output and input
using namespace std;

// object for output files
ofstream ofile;
// Functions used
inline double f(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

// Begin main program
int main(int argc, char *argv[]){
  int exponent; 
    string filename;
    // We read also the basic name for the output file and the highest power of 10^n we want
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max power 10^n" << endl;
          exit(1);
    }
        else{
        filename = argv[1]; // first command line argument after name of program
        exponent = atoi(argv[2]);
    }
    // Loop over powers of 10
    for (int i = 1; i <= exponent; i++){
      int  n = (int) pow(10.0,i);
      // Declare new file name
      string fileout = filename;
      // Convert the power 10^i to a string
      string argument = to_string(i);
      // Final filename as filename-i-
      fileout.append(argument);
      double h = 1.0/(n);
      double hh = h*h;
      // Set up arrays for the simple case
      double *d = new double [n+1]; double *b = new double [n+1]; double *solution = new double [n+1];
      double *x = new double[n+1];
      // Quick setup of updated diagonal elements and value of
      d[0] = d[n] = 2; solution[0] = solution[n] = 0.0;
      for (int i = 1; i < n; i++) d[i] = (i+1.0)/( (double) i);  
      for (int i = 0; i <= n; i++){
	x[i]= i*h;
        b[i] = hh*f(i*h);
      }
      // Forward substitution
      for (int i = 2; i < n; i++) b[i] = b[i] + b[i-1]/d[i-1];
      // Backward substitution
      solution[n-1] = b[n-1]/d[n-1];
      for (int i = n-2; i > 0; i--) solution[i] = (b[i]+solution[i+1])/d[i];
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i < n;i++) {
	double xval = x[i];
 	 double RelativeError = fabs((exact(xval)-solution[i])/exact(xval));
         ofile << setw(15) << setprecision(8) << xval;
         ofile << setw(15) << setprecision(8) << solution[i];
         ofile << setw(15) << setprecision(8) << exact(xval);
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();
      delete [] x; delete [] d; delete [] b; delete [] solution;
    }
    return 0;
}


