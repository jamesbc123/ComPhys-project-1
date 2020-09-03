#include <iostream>
#include "gaussElim.h"

using namespace std

double x0 = 0;  // Starting point
double x_np1 = 1;  // End point (x_(n+1))
// Boundary conditions:
double u0 = 0; // u_0 = u(x0) = u(0)
double u_np1 = 0; // u_(n+1) = u(xn) = u(1)

int n = 1000;  // Number of unknowns (u1, u2, ... , un).

//test comment
