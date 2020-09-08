#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include "gaussElim.cpp"

using namespace std;

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
    int p = 7; // Power of 10.
    int n = pow(10, p);
    cout << "n: " << n << endl;

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
    for (int i=0; i<=n-2; i++) {
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
    }
    b[n-1] = 2; // Fill the last element of b.


    double* solutionGen = new double[n+2];
    double* solutionSpec = new double[n+2];

    double tGen = general_algo(solutionGen, n, x_0, x_np1, a, b, c);  // Time elapsed for the 
    // general algorithm.
    double tSpec = special_algo(solutionSpec, n, x_0, x_np1, a[0], b[0]); // Time elapsed for the
    // special algorithm.

    cout << "tGen: " << tGen << ". tSpec: " << tSpec << endl;

    // ### STEP 3 ### Choose what values to write to file (what to plot).
    // Choice: Plot the solution as a function of x:
    /*double* argumentList = new double[n+2];
    double* valueList = new double[n+2];
    argumentList = xList;   // The x-axis of the plot.
    valueList = solution;   // The y-axis of the plot. */

    //printArray(xList, n+2);
    //printArray(solution, n+2);

    // ### STEP 4 ### Write the solution to a .csv file.
    string fileNameGen = "data_general.txt";
    string fileNameSpec = "data_special.txt";
    int listLength = n+2;

    // Time writing to file.
    writeDataToCSV(xList, solutionGen, listLength, fileNameGen);
    writeDataToCSV(xList, solutionSpec, listLength, fileNameSpec);

    // Delete all allocated memory:
    delete[] a; delete[] b; delete[] c;
    delete[] xList; delete[] solutionGen; delete[] solutionSpec;
    return 0;

    // ### STEP 5 ### In Python: Read the .csv file and plot it.
}
