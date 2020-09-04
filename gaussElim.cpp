#include <iostream>
#include <cmath>
#include "time.h"
#include "gaussElim.h"

double* general_algo(int n, double x_0, double x_np1, , a, b, c)
{
    /* 
    Returns solution to general algo and time elapsed.

    This function calculates v for the general case, where
    a, b and c are different values. This is done by forward 
    and backward substitution. 

    Inputs: 
        a, b, c are vectors.
        a: Lower diagonal elements (n-1 elements).
        b: Middle diagonal (n elements).
        c: Upper diagonal (n-1 elements).
        n: the number of points excluding the two end points (in total n+2 points).
        x_0: Start point.
        x_np1: End point (x_(n+1)).
        fx: The values of f(x) for 
    Outputs:
        v: solution, array of length n+2.
    */

    // Calculate h. Using long double for increased precision. (Does *everything* have to be long double)
    // in order to get less numerical errors?)
    long double h = (x_np1-x_0)/(n+1);
    long double hh = h*h;

    // Create xList from end points (x_0, x_np1).
    idx = 0;
    double *xList = new double[n+2]; // n+2 points in total when including the end points. 
    xList[0] = x_0; // End points.
    xList[n+1] = x_np1;
    for (i=1; i<=n; i++) { xList[i] = x_0 + i*h; }    // Step size h between points.

    // Create fList (an array containing f(x) for all elements in xList, including end points):
    double *fList = new double[n+2];
    double *gList = new double[n+2];

    double *solution = new double[n+2]; // Solution for v = [v0, v1, v2, ..., v_n, v_(n+1)].
    solution[0] = 0; solution[n+2] = 0; // Boundary conditions.

    fList[0] = f(x_0); // End points (both f and g arrays).
    gList[0] = hh*fList[0];

    fList[n+1] = f(x_np1);
    gList[n+1] = hh*fList[n+1];     // Remember: g_i = h^2*f_i  (g_i is b_tilde_i in the project 1 description pdf.)

    for (i=0; i<=n+1; i++) {
        fList[i] = f(xList[i]);
        gList[i] = hh*fList[i];
    }

    // Now: The set of linear equations we need to solve is Av = g (solve for the vector v), for a tridiagonal matrix A.
    // The set of equations will not be solved by creating the actual matrix and going through each row, but just by
    // taking the diagonal elements as inputs and creating each matrix row 'on the go' inside the for-loop. (Because 
    // creating the nxn matrix is computationally expensive.) 
    // It is solved with gaussian elimination, which is first forward substitution and then backward substitution.

    // Forward substitution (start at row 2 and to all subsequent rows).
    
    for (i=1; i<=n-1; i++) { // Forward substitution goes from row 2 (index 1) through row n (index n-1).
        // The matrix A is defined by the diagonal vectors a, b, c.
        // Note: the lower diagonal, a, starts at row 2, and the middle and upper diagonals, b and c, starts at row 1.
        // This means that a[i] is on row i and b[i] and c[i] are on row i-1. 

        // Pre-calculate factor:
        double factorDiv = b[i-1]; // Divide the row by this number before subtracting (Gaussian elimination).
        // (The element directly above the leading element of the current row.)
        double factorMult = a[i]; // Multiply the row by this number before subtracting. (Leading element of the current row.)
        double factor = factorMult/factorDiv; // This is the factor to be multiplied to the row before subtracting.
=
        // First non-zero column term (on current row, row i):
        //a[i-1] = a[i-1] - b[i-1]*f; // We don't need 'a' for any later calculations, so this calculation is omitted.
        
        // Second non-zero column term:
        b[i] = b[i] - c[i-1]*f; // 'b' is used in the backward substitution, so this calculation is necessary.
        
        // Since the matrix is tri-diagonal, the third non-zero column element has a zero directly above it,
        // so nothing is subtracted. So nothing happens here. c[i] = c[i] - 0.

        // Final column element:
        //M[i,-1] = M[i,-1] - M[i-1,-1]*f
        g[i] = g[i] - g[i-1]*factor;
    }

    // Backward substitution and normalization (and the solution):
    for (i=n-2; i>=0; i--) { // Backward substitution goes from row n-1 (index n-2) and row 1 (index 0).
        // Pre-calculate factor:
        double factorDiv = b[i+1];
        double factorMult = c[i];
        double factor = factorMult/factorDiv;

        // First non-zero column element is eliminated:
        //c[i] = c[i] - b[i+1]*factor =0 

        // Final column element:
        g[i] = g[i] - g[i+1]*factor;
        // Normalize the diagonal (b) in order to get the solution:
        g[i] = g[i]/b[i];
        solution[i] = g[i];
    }

    // Return solution:
    return solution;
}


double* special_algo(int n, double x_0, double x_np1, a, b)
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
delete [] b_sub; delete [] g_sub; delete [] v;

return v, ((finish - start)/CLOCKS_PER_SEC); 
// End special algo function
}