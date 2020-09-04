import numpy as np
import matplotlib.pyplot as plt

# This is the code for project 1 in FYS4150, computational physics.

def createTriDiagonalMatrix(N, a, b, c):
    # This function returns a tri-diagonal matrix with specified values.
    # N is the dimension. The output matrix will be an NxN matrix. N should be equal to or greater than 3.
    # Elements on the below-diagonal will have the value a,
    # elements on the diagonal will have the value b,
    # elements on the above-diagonal will have the value c.
    M = np.zeros((N,N))

    # Note: Python uses 0-indexing.

    # The first row:
    M[0,0] = b
    M[0,1] = c
    # The final row:
    M[N-1,N-2] = a
    M[N-1,N-1] = b

    # The other rows:
    for i in range(1,N-1):
        M[i, i-1] = a
        M[i, i] = b
        M[i, i+1] = c
    return M

# IMPORTANT: Since M is a tridiagonal matrix it's not necessary to create the matrix
# beforehand. Creating the matrix is very memory costly since its nxn (n^2). It's better to
# just take each diagonal as inputs (as vectors, ~n elements each). The C++ translation
# of this program will do it this way.
def createAugmentedMatrix(M, vec):
    # This appends the vector vec as the final column of M.
    # Note: The length of vec must be compatible with M. In other words, the length
    # of vec must be equal to the number of rows of M.
    dim = np.array(np.shape(M)) # Matrix dimensions.

    newDimensions = dim + [0,1] # Add one column to M.
    augM = np.zeros(newDimensions)
    augM[:,:-1] = M
    augM[:,-1] = vec    # Set final column to vec.
    return augM


def gaussForwardSubstTriDiag(M):
    # This function performs the first part of gaussian elimination, called forward substitution, on the
    # system of equations M, equivalent to Ax = b, where A is a tri-diagonal matrix.
    # M has dimensions N x (N+1).

    N = np.shape(M)[0]  # The number of rows of M.
    for i in range(1,N):   # Go through and subtract the rows, starting with the second row.
        f1 = M[i-1,i-1]    # Divide the row by this number before subtracting (Gaussian elimination).
        # (The element directly above the leading element of the current row.)
        f2 = M[i,i-1]    # Multiply the row by this number before subtracting. (Leading element of the current row.)
        f = f2/f1   # This is the factor to be multiplied to the row before subtracting.

        # Subtract (a multiple of) the previous row (i-1) from the current row (i) in order to eliminate
        # the leading term in the current row:
        # First non-zero column term:
        M[i,i-1] = M[i,i-1] - M[i-1,i-1]*f
        # Second non-zero column term:
        M[i,i] = M[i,i] - M[i-1,i]*f
        # Since the matrix is tri-diagonal the third non-zero column element of the current row has a zero in the
        # element in the row above it, so nothing is subtracted. So M[i,i+1] remains the same.
        # Final column element:
        M[i,-1] = M[i,-1] - M[i-1,-1]*f

    return M


def gaussBackwardSubstTriDiag(M):
    # This function performs the second part of gaussian elimination, called backward substitution.
    # The input matrix M has to be an upper triangular matrix with non-zero elements only along the diagonal
    # and on the adjacent diagonal above it. This function simply eliminates the terms on this other diagonal.

    N = np.shape(M)[0]  # The number of rows of M.
    for i in range(N-2,-1,-1):
        f1 = M[i+1,i+1]    # Divide the row by this number before subtracting (Gaussian elimination).
        # (The element directly above the leading element of the current row.)
        f2 = M[i,i+1]    # Multiply the row by this number before subtracting. (Leading element of the current row.)
        f = f2/f1   # This is the factor to be multiplied to the row before subtracting.

        # Eliminate the upper diagonal element:
        M[i,i+1] = M[i,i+1] - M[i+1,i+1]*f

        # Final column element:
        M[i,-1] = M[i,-1] - M[i+1,-1]*f

    return M

def triDiagonalgaussElim(A, bVec):
    # This function returns the vector x that is the solution of the system of equations Ax = b,
    # which is represented by the augmentet matrix M. A is a tridiagonal matrix.
    # The system is solved by Gaussian elimination.

    # Create the augmented matrix for the system of linear equations:
    M = createAugmentedMatrix(A, bVec)

    # Forward substitution:
    M = gaussForwardSubstTriDiag(M)
    # Backward substitution:
    M = gaussBackwardSubstTriDiag(M)

    # Normalize the diagonal so that we get the solution:
    N = np.shape(M)[0]  # The number of rows of M.
    for i in range(0,N):    # Go through all the rows.
        # Final column:
        M[i,-1] = M[i,-1]/M[i,i]
        # Diagonal elements:
        M[i,i] = 1; # In reality it says M[i,i] = M[i,i]/M[i,i], but we know that this is 1.

    # Extract the solution vector (the final column):
    xSolution = M[:,-1]
    return xSolution

def f(x):
    # This is the 'source term' in Poissons equation -u''(x) = f(x) in this case.
    # Physically f(x) would correspond to a spherically symmetric charge distribution,
    # e.g. a point charge, while the unknown u(x) would be the resulting electric
    # potential caused by this charge distribution, and x is the distance from its
    # center.
    return 100*np.e**(-10*x)

def uAnalytical(x):
    # This is the analytical solution for u of the equation -u''(x) = f(x) in the
    # case that f(x) = 100e^(-10x).
    return 1 - (1-np.e**(-10))*x - np.e**(-10*x)

def plotNumericalVsAnalytical(xList, numericalSol):
    # Plots the numerical solution vs the analytical solution (u(x)) of the Poisson equation.
    n = len(xList) - 2  # Remove the two end points.
    xListA = np.linspace(x0, x_np1, 1000)   # Plot the analytical solution with high resolution.
    analyticalSolHR = uAnalytical(xListA)  # The exact analytical solution for u(x) (high-resolution)
    # the numerical solution (in order to compare them).

    plot1, = plt.plot(xListA, analyticalSolHR, label = 'Analytical', linewidth=2)
    plot2, = plt.plot(xList, numericalSol, '-r', label = 'Numerical', linewidth=1)

    plt.legend()
    plt.grid()
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    #plt.xlim(0,max(xliste)*1.05)  # Setter grensene for x-aksen.
    #plt.ylim(0,max(yliste)*(1.05)) # Setter grensene for y-aksen.
    plt.suptitle("Analytical vs. numerical solution for " + r'$u(x)$' + " with n = %d" %n)
    plt.show()

def plotErrors(xList, numericalSol):
    # Plots the logarithms of the absolute values of the relative errors for all
    # calculated points.
    # numericalSol: The numerical solution of u(x) for the values
    # of x in xList.
    numericalSol = numericalSol[1:-1] # Remove the end points.
    xList = xList[1:-1] # Remove the end points.
    n = len(xList)
    analyticalSol = uAnalytical(xList)
    print(numericalSol)
    print(analyticalSol)
    print()

    errors = numericalSol - analyticalSol
    errors = abs(errors)
    errorsD = errors/analyticalSol
    print(errors)
    print(errorsD)
    print()
    logErrors = np.log10(abs(errors/analyticalSol)) # Element-wise log_10 of the
    # absolute values of the relative errors.
    print(logErrors)

    #plt.plot(xList, relativeErrors, label = r'$$\log_{10} (|\frac{v_i-u_i}{u_i}|)$$', linewidth=2)
    #plt.plot(xList, relativeErrors, label = 'error', linewidth=2)

    #plt.plot(xList, errorsD, label = 'error')
    #plt.yscale('log')
    #plt.semilogy(xList, errorsD)#, label = 'error', linewidth=1)

    #plt.plot(xList, errors, label = r'Error: $\left|v_i-u_i\right|$')

    plt.plot(xList, logErrors, label = r'Error: $\epsilon_i = \log\left(\left|\frac{v_i-u_i}{u_i}\right|\right)$')
    plt.gca().set_ylim(-6, 0)

    plt.legend()
    plt.grid()
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.suptitle("Numerical error for n = %d" %n)
    plt.show()



# ##############################

x0 = 0  # Starting point
x_np1 = 1  # End point (x_(n+1))
# Boundary conditions:
u0 = 0 # u_0 = u(x0) = u(0)
u_np1 = 0 # u_(n+1) = u(xn) = u(1)

n = 10  # Number of unknowns (u1, u2, ... , un).

# The step length is defined by the interval and the number of grid points:
# There are n+2 grid points (since the edges, x0 and x_np1, are included).
h = (x_np1-x0)/(n+1)   # With n+2 grid points there are n+1 sub-intervals.
print(h)
print()

A = createTriDiagonalMatrix(n, -1, 2, -1)
#print(A)
#print()

# Now get the numerical solution:
xList = np.linspace(x0, x_np1, n+2)    # n+2 grid points
bTildeVec = h**2 * f(xList[1:-1]) # The b-tilde vector from the equation Av = b-tilde.
# ^ The first and last elements in xList are not included in the set of linear equations.
numericalSol = np.zeros(n+2)
numericalSol[0], numericalSol[-1] = 0, 0    # Boundary conditions.
numericalSol[1:-1] = triDiagonalgaussElim(A, bTildeVec)

#plotNumericalVsAnalytical(xList, numericalSol)
plotErrors(xList, numericalSol)

