#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include "gaussElim.h"

/*
double* general_algo(int n, double x_0, double x_np1, double* a, double* b, double* c)
{

    Returns solution to general algo and time elapsed.

    This function calculates v for the general case, where
    a, b and c are different values. This is done by forward 
    and backward substitution. 
    
    Inputs: 
        a, b, c are arrays containing the elements of the diagonals.
        a has length n, b has length n+1, c has length n.
        n: the number of points excluding the two end points (in total n+2 points).
        x_0: Start point.
        x_np1: End point (x_(n+1)).
    Outputs:
        v: solution, array of length n+2.


    double h = (x_np1-x_0)/(n+1);
    double hh = h*h;

    // Create xList from end points (x_0, x_np1).
    double *xList = new double[n+2]; // n+2 points in total when including the end points. 
    xList[0] = x_0; // End points.
    xList[n+1] = x_np1;
    for (int i=1; i<=n; i++) { xList[i] = x_0 + i*h; }    // Step size h between points.

    // Create fList (an array containing f(x) for all elements in xList, including end points):
    double *fList = new double[n+2];
    double *g = new double[n+2];    // g_i = h^2*f_i
    double *solution = new double[n+2]; // The solution of v.
    solution[0] = 0; solution[n+1] = 0; // Boundary conditions.

    for (int i=0; i<=n+1; i++) { // Including end points.
        fList[i] = f(xList[i]);
        g[i] = hh*fList[i];
    }

    // Forward substitution:
    for (int i=1; i<=n-1; i++) { // Start at row 2 (index 1), end at final row (row n, index n-1).
        double factorDiv = b[i-1];
        double factorMult = a[i-1];
        double factor = factorMult/factorDiv;

        // First non-zero element in the row:
        b[i] = b[i] - c[i-1]*factor;
        // Final element in the row:
        g[i] = g[i] - g[i-1]*factor;
    }

    // Backward substitution:
    for (int i=n-2; i>=0; i--) { // Start at row n-1 (index n-2), end at first row (row 1, index 0).
        double factorDiv = b[i+1];
        double factorMult = c[i];
        double factor = factorMult/factorDiv;
        // All upper diagonal elements gets eliminated.
        // Final element in the row:
        g[i] = g[i] - g[i+1]*factor;
    }

    // Normalize the diagonal (divide all row i by b[i] for all rows) in order to get the 
    // solution for v:
    for (int i=0; i<=n-1; i++) {
        solution[i] = g[i]/b[i];
    }
    delete[] fList; delete[] g; delete[] xList; 
    return solution;
}

*/

void general_algo(double* solution, int n, double x_0, double x_np1, double* a, double* b, double* c)
{
    /* 
    Returns solution to general algo and time elapsed.

    This function calculates v for the general case, where
    a, b and c are different values. This is done by forward 
    and backward substitution. 
    
    Inputs: 
        solution: A dynamic array of length n+2 allocated outside of this function.
        This function fills this array with the solution values.
        a, b, c are arrays containing the elements of the diagonals.
        a has length n, b has length n+1, c has length n.
        n: the number of points excluding the two end points (in total n+2 points).
        x_0: Start point.
        x_np1: End point (x_(n+1)).
    Outputs:
        No direct outputs, but the values at the address of the 'solution' array will
        be filled with the solution values.
    */

    double h = (x_np1-x_0)/(n+1);
    double hh = h*h;

    // Create xList from end points (x_0, x_np1).
    double *xList = new double[n+2]; // n+2 points in total when including the end points. 
    xList[0] = x_0; // End points.
    xList[n+1] = x_np1;
    for (int i=1; i<=n; i++) { xList[i] = x_0 + i*h; }    // Step size h between points.

    // Create fList (an array containing f(x) for all elements in xList, including end points):
    double *fList = new double[n+2];
    double *g = new double[n+2];    // g_i = h^2*f_i
    //double *solution = new double[n+2]; // The solution of v.
    solution[0] = 0; solution[n+1] = 0; // Boundary conditions.

    for (int i=0; i<=n+1; i++) { // Including end points.
        fList[i] = f(xList[i]);
        g[i] = hh*fList[i];
    }

    // Forward substitution:
    for (int i=1; i<=n-1; i++) { // Start at row 2 (index 1), end at final row (row n, index n-1).
        double factorDiv = b[i-1];
        double factorMult = a[i-1];
        double factor = factorMult/factorDiv;

        // First non-zero element in the row:
        b[i] = b[i] - c[i-1]*factor;
        // Final element in the row:
        g[i] = g[i] - g[i-1]*factor;
    }

    // Backward substitution:
    for (int i=n-2; i>=0; i--) { // Start at row n-1 (index n-2), end at first row (row 1, index 0).
        double factorDiv = b[i+1];
        double factorMult = c[i];
        double factor = factorMult/factorDiv;
        // All upper diagonal elements gets eliminated.
        // Final element in the row:
        g[i] = g[i] - g[i+1]*factor;
    }

    // Normalize the diagonal (divide all row i by b[i] for all rows) in order to get the 
    // solution for v:
    for (int i=0; i<=n-1; i++) {
        solution[i] = g[i]/b[i];
    }
    delete[] fList; delete[] g; delete[] xList; 
    //return solution;
}




double* special_algo(int n, double x_0, double x_np1, int a, int b) {
    /* Returns solution to the special algorithm and elapsed time.

    This function calculates v for the general case, i.e a 
    tri-diagonal matrix, where as opposed to the general case
    a and c are equal and all elements along the diagonal and 
    off-diagonals are identical. 
    Inputs:
        a, b: integers, diagonal elements
        n: interger, number of points
        x_0: double, point 0
        x_np1: double, point n+1
    Outputs:
        v: double, solution of size n+1 */

    // Calculate h. Using long double for increased precision.
    double h = (x_np1-x_0)/(n+1);
    double hh = h*h;
    double *fList = new double[n+2]; double *g = new double[n+2]; 
    double *xList = new double[n+2]; double *v = new double[n+2];
    double *b_sub = new double[n+2]; double *g_sub = new double[n+2];

    // Create xList:
    xList[0] = x_0; // End points.
    xList[n+1] = x_np1;

    fList[0] = x_0; // End points.
    fList[n+1] = x_np1;
    for (int i=1; i<=n; i++) { xList[i] = x_0 + i*h; }    // Step size h between points.

    for (int i=0; i<=n+1; i++) {
        fList[i] = f(xList[i]);
        g[i] = hh * fList[i];
    }

    //clock_t start, finish; // declare start and final time
    //start = clock();

    // Forward substitution
    b_sub[1] = b;
    for (int i=2; i<=n+1; i++){
        b_sub[i] = (i+1)/i;
        g_sub[i] = g[i] + (g_sub[i-1]/b_sub[i-1]);
        } 
    v[n] = g_sub[n] / b_sub[n];

    // Backward substitution
    for (int i=n-1; i>0; i--){
        v[i] = (g_sub[i] + v[i+1])/b_sub[i];
        }
            
    //finish = clock();

    delete [] xList; delete [] fList; delete [] g;
    delete [] b_sub; delete [] g_sub;

    return v;//, ((finish - start)/CLOCKS_PER_SEC);
}

double* relative_error(double* v, double* u, int n){
    /* Returns the relative error between v and u.
    
    This function compares the exact value to the numerical solution.
    Inputs:
        u: an array of the exact values.
        v: an array of numerical values.
        n: Number of points being compared.
    Outputs:
        rel_err: an array of relative errors.
    */
    
    double *rel_err = new double [n+2];
    for (int i=1; i<=n; i++){
        rel_err[i] = log10(fabs((u[i] - v[i])/u[i]));
    }
    return rel_err;
}


void writeDataToCSV(double* xList, double* yList, int listLength, std::string fileName) {
    /* Writes the info in xList and yList as two separate columns to a comma-separated .txt
    file. This way the data can be plotted using e.g. Python.
    xList: Values along the x-axis.
    yList: The corresponding data values along the y-axis. Same length as xList.
    listLength: The number of elements (length) of xList and yList.
    file_name: The chosen name of the new data file. If this already exists, the old file
    will be overwritten.
    */
    std::ofstream file;
    file.open(fileName);
    // file.open(fileName, std::ios_base::app); // Use this instead of you want the program
    // to *append* to the file instead of overwriting.

    for (int i=0; i<=listLength-1; i++) {
        file << std::to_string(xList[i]) << "," << std::to_string(yList[i]) << "," << std::endl;
    }
    
    file.close();
}

