#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include "gaussElim.h"

using namespace std;
// object for output files
ofstream ofile;


/* Description of the whole data plotting process:
- Set the initial parameter values of the problem.
- Use the functions in gaussElim.h to solve the system of linear equations (using either
general_algo or special_algo).
- Write this data to a .csv file using the writeDataToCSV function.
- Open the Python file main.py and read the data from the .csv file.
- Plot the data using Python. Make sure the title and labels are correct
for each plot.
*/
int main() {
    // This case: General algorithm.

    // ### STEP 1 ### Initialize the parameters of the problem.
    // Boundaries for x:
    double x_0 = 0;
    double x_np1 = 1;

    // Choose the number of grid points, n:
    int n = 10;

    // Step value:
    double h = (x_np1-x_0)/(n+1);

    // Create xList, the list of all x values (n+2 points including end points):
    double *xList = new double[n+2];
    xList[0] = x_0;
    xList[n+1] = x_np1;
    for (int i=1; i<=n; i++) { xList[i] = x_0 + i*h; } 

    // ### STEP 2 ### Get the solution of the differential equation.
    // Make the diagonal arrays:
    double *a = new double[n-1];    // Lower diagonal.
    double *b = new double[n];      // Middle diagonal.
    double *c = new double[n-1];    // Upper diagonal.
    for (int i=0; i<=n-1; i++) {
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
    }

    double* solution = general_algo(n, x_0, x_np1, a, b, c);
    // ^ 'solution' is an array of length n+2.

    // ### STEP 3 ### Choose what values to write to file (what to plot).
    // Choice: Plot the solution as a function of x:
    double* argumentList = xList;   // The x-axis of the plot.
    double* valueList = solution;   // The y-axis of the plot.

    // ### STEP 4 ### Write the solution to a .csv file.
    string fileName = "data.txt"; // Or 'data.csv'?
    writeDataToCSV(argumentList, valueList, fileName);

    // Delete all allocated memory:
    delete[] xList; delete[] a; delete[] b; delete[] c;
    delete[] argumentList; delete[] valueList;

    return 0;

    // ### STEP 5 ### In Python: Read the .csv file and plot it.
}


/*
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

    double x_0 = 0;  // Starting point
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
      double *v = new double[n+2];
      double time_special;
      v_special, time_special = special_algo(n, x_0, x_np1, a , b);
      double *x = new double[n+2];
      double *u = new double[n+2];
      x[0] = x_0; // End points.
      x[n+1] = x_np1;
      
      for (i=1; i<=n; i++){
          x[i] = x_0 + i*h;
          u[i] = exact(x[i]);
      }
      double *rel_err = new double [n+2];
      rel_err = relative_error(v_special, u);

      ofile << setw(15) << setprecision(8) << rel_err << endl;
      
      ofile.close();
      
    }
    return 0;
}
*/