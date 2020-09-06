#ifndef GAUSSELIM_H
#define GAUSSELIM_H

#include <iostream>
#include <cmath>

// Declare all the functions to be used for Gaussian elimination here,
// and implement those functions in gaussElim.cpp


double f(double x) {return 100.0*exp(-10.0*x);}

double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double* general_algo(int n, double x_0, double x_np1, double* a, double* b, double* c);
/* Returns the solution to the general algorithm.
This function calculates v for the general case, where the diagonals
a, b and c can have any elements. This is done by forward and backward
substitution. 
Inputs: 
    a, b, c are arrays containing the elements of the diagonals.
    a has length n, b has length n+1, c has length n.
    n: the number of points, excluding the two end points (in total: n+2 points).
    x_0: Start point.
    x_np1: End point (x_(n+1)).
Outputs:
    v: solution, array of length n+2.
*/


double* special_algo(int n, double x_0, double x_np1, int a, int b);
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
    v: double, solution of size n+1
*/

double* relative_error(double v, double u);
/* Returns the relative error between v and u.

This function compares the exact value to the numerical solution.
Inputs:
    u: an array of the exact values.
    v: an array of numerical values.
Outputs:
    rel_err: an array of relative errors.
*/


void writeDataToCSV(double* xList, double* yList, std::string fileName);
/* Writes the info in xList and yList as two separate columns to a comma-separated .txt
file. This way the data can be plotted using e.g. Python. 
xList: Values along the x-axis.
yList: The corresponding data values along the y-axis.
file_name: The chosen name of the new data file. If this already exists, the old file
will be overwritten.
*/
    



#endif // GAUSSELIM_H