#ifndef GAUSSELIM_H
#define GAUSSELIM_H

#include <iostream>
#include <cmath>

// declare all the functions to be used for Gauss elimination here,
// and implement those functions in gaussElim.cpp


double f(double x) {return 100.0*exp(-10.0*x);}

double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double* general_algo(int n, double x_0, double x_np1, double a, double b, double c);
/* Returns solution to general algo and time elapsed.
This function calculates v for the general case, where
a, b and c are different values. This is done by forward 
and backward substitution. 
Inputs: 
    a, b, c are arrays of length n+1
    n: the nbr of points
    hh: the step length squared
Outputs:
    v: solution, array of length n+2 */


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
    v: double, solution of size n+1 */

double* relative_error(double v, double u);
    /* Returns the relative error between v and u.
    
    This function compares the exact value to the numerical solution.
    Inputs:
        u: an array of the exact values.
        v: an array of numerical values.
    Outputs:
        rel_err: an array of relative errors.
    */




#endif // GAUSSELIM_H