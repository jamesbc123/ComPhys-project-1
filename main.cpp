#include <iostream>
#include "gaussElim.h"

using namespace std


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

    double x0 = 0;  // Starting point
    double x_np1 = 1;  // End point (x_(n+1))
    
    // Boundary conditions:
    double u0 = 0; // u_0 = u(x0) = u(0)
    double u_np1 = 0; // u_(n+1) = u(xn) = u(1)
    int a = -1;
    int b = 2;

    // Loop over powers of 10
    for (int i = 1; i <= exponent; i++){
      int  n = (int) pow(10.0,i);
      double h = (x_np1-x_0)/(n+1);

      // Declare new file name
      string fileout = filename;

      // Convert the power 10^i to a string
      string argument = to_string(i);

      // Final filename as filename-i-
      fileout.append(argument);
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      ofile << "       x:             approx:          exact:       relative error" << endl;

      vSpecial, timeSpecial = special_algo(n, x_0, x_np1, a , b);
      double *x = new double[n+2];
      double *u = new double[n+2];
      x[0] = x_0; // End points.
      x[n+1] = x_np1;
      
      for (i=1; i<=n; i++){
          x[i] = x_0 + i*h;
          u[i] = exact(x[i]);
      }
      double rel_err = relative_error(vSpecial, u);

      ofile << setw(15) << setprecision(8) << rel_err << endl;
      
      ofile.close();
      
    }
    return 0;
}
