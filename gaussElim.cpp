#include <iostream>
#include <cmath>
#include "gaussElim.h"


double* general_algo(int n, double x_0, double x_np1, , a, b, c)
{
    /* 
    Returns solution to general algo and time elapsed.

    This function calculates v for the general case, where
    a, b and c are different values. This is done by forward 
    and backward substitution. 

    Inputs: 
        a, b, c are arrays of length n+1
        n: the number of points excluding the two end points (in total n+2 points).
        x_0: Start point.
        x_np1: End point (x_(n+1)).
        fx: The values of f(x) for 
    Outputs:
        v: solution, array of length n+2.
    */
   x
    // Calculate h. Using long double for increased precision.
    long double h = (x_np1-x_0)/(n+1);

    // Create xList from end points (x_0, x_np1).
    idx = 0;
    double *xList = new double[n+2]; // n+2 points in total when including the end points. 
    xList[0] = x_0; // End points.
    xList[n+1] = x_np1;
    for (i=1; i<=n; i++) { xList[i] = x_0 + i*h; }    // Step size h between points.

    // Create fList (an array containing f(x) for all elements in xList, including end points):
    double *fList = new double[n+2];
    fList[0] = x_0; // End points.
    xList[n+1] = x_np1;
    for (i=0; i<=n+1; i++) {
        fList[i] = f(xList[i]);
    }

    // FERDIG MED fList, FORTSETT MED FORWARD SUBSTITUTION ******************************

    // Forward substitution:


    // Backward substitution:


    // Return solution:

}


double* special_algo(int n, double x_0, double x_np1, , a, b, c)
{

}