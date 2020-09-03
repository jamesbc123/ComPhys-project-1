#ifndef GAUSSELIM_H
#define GAUSSELIM_H

#include <iostream>
#include <cmath>

// declare all the functions to be used for Gauss elimination here,
// and implement those functions in gaussElim.cpp


double f(double x) {return 100.0*exp(-10.0*x);}

double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double* general_algo(int n, double x_0, double x_np1, , a, b, c);
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


double* special_algo(n, h);






#endif // GAUSSELIM_H