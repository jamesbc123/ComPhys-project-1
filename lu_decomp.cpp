#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include "gaussElim.cpp"

using namespace std;   
using namespace arma;
ofstream ofile;

int main()
{

int n = 10;
double x_0 = 0;
double x_np1 = 1;
double h = (x_np1-x_0)/(n+1);
double hh = h*h;
double *x = new double[n];
mat u = zeros<vec>(n);
mat fx = zeros<vec>(n);
mat solution;
x[0] = h;

for(int i=0; i<n; i++) {
    x[i] = x[0] + i*h;
    fx[i] = f(x[i]) * hh;
    u[i] = exact(x[i]);
    }

// Calculate the solution using armadillo package
mat A = zeros<mat>(n, n);
A.diag().fill(2);
A.diag(-1).fill(-1);
A.diag(1).fill(-1);

// Lower matrix, Upper matrix and Permutation matrix.
mat L, U, P;
lu(L, U, P, A);

// Provide a permutation such that P.t()*L*U = A.
//mat B = P.t()*L*U;

// To solve Ax=(LU)x=b first solve Ly=b for y
mat y = solve(L, fx);

// now solve for x
solution = solve(U, y);

double* relative_err = new double[n];
for(int i=0; i<n; i++){
        relative_err[i] = log10(abs((u[i] - solution[i])/u[i]));
}

string fileout = "lu_decomp.txt";
ofile.open(fileout);

for(int i=0; i<n; i++){
ofile << setw(15) << setprecision(8) << x[i];
ofile << setw(15) << setprecision(8) << solution[i];
ofile << setw(15) << setprecision(8) << u[i];
ofile << setw(15) << setprecision(8) << relative_err[i] << endl;
}
ofile.close();

return 0;
}

