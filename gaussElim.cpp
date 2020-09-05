#include <iostream>
#include <cmath>
#include "gaussElim.h"
#include "time.h"

double* general_algo(int n, double x_0, double x_np1, double a, double b, double c)
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


double* special_algo(int n, double x_0, double x_np1, int a, int b)
{
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
double hh = h*h
double *fList = new double[n+2]; double *g = new double[n+2]; 
double *xList = new double[n+2]; double *v = new double[n+2];
double *b_sub = new double[n+2]; double *g_sub = new double[n+2];

// Create xList from end points (x_0, x_np1).
idx = 0;
xList[0] = x_0; // End points.
xList[n+1] = x_np1;

fList[0] = x_0; // End points.
fList[n+1] = x_np1;
for (i=1; i<=n; i++) { xList[i] = x_0 + i*h; }    // Step size h between points.

for (i=0; i<=n+1; i++) {
    fList[i] = f(xList[i]);
    g[i] = hh * fList[i];
}

clock_t start, finish; // declare start and final time
start = clock();

// Forward substitution
b_sub[1] = b 
for (i=2; i<=n+1; i++){
    b_sub[i] = (i+1)/i
    g_sub[i] = g[i] + (g_sub[i-1]/b_sub[i-1])
    } 
v[n] = g_sub[n] / b_sub[n]

// Backward substitution
for (i=n-1; i>0; i--){
    v[i] = (g_sub[i] + v[i+1])/b_sub[i]
    }
        
finish = clock();

delete [] xList; delete [] fList; delete [] g;
delete [] b_sub; delete [] g_sub;

return v, ((finish - start)/CLOCKS_PER_SEC); 
// End special algo function
}

double* relative_error(double v, double u){
    /* Returns the relative error between v and u.
    
    This function compares the exact value to the numerical solution.
    Inputs:
        u: an array of the exact values.
        v: an array of numerical values.
    Outputs:
        rel_err: an array of relative errors.
    */
    
    double *rel_err = new double [n+2];
    for (i=1, i<=n, i++){
        rel_err[i] = log10(fabs((u[i] - v[i])/u[i]))
    }
    return rel_err;
}