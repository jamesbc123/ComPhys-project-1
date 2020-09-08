#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include "gaussElim.cpp"

using namespace std;

/*
Program description:
This program times the general algorithm and the special algorithm.
The average run times for different values of n are written to .csv
files, which can then be plotted using Python.
*/

int main() {
    // Iterate over powers of 10: 10, 100, ... , 10^6, 10^7.
    int powerMax = 7;   // Maximum power of 10.
    int listLength = powerMax;  // Assuming the powers are 1, 2, ..., powerMax.
    double* nList = new double[listLength];  // x axis in plot (n).
    double* timeListGeneral = new double[listLength]; // y axis in plot (time).
    double* timeListSpecial = new double[listLength];

    int numAverage = 4;    // Record and calculate the *average* CPU time over this many
    // runs (for each n).

    for (int p=1; p<=powerMax; p++) {
        cout << "Power of 10: n = " << p << "." << endl;
        int n = pow(10, p); // n = 10^p
        nList[p-1] = n;
        
        // Total time for the general and special case (sum of all iterations).
        double tSumGen = 0;
        double tSumSpec = 0;

        for (int iter = 1; iter <= numAverage; iter++) { // Time the general and special algorithms numAverage times.
            double* solutionGeneral = new double[n+2];
            double* solutionSpecial = new double[n+2];

            // ### STEP 1 ### Initialize the parameters of the problem.
            // Boundaries for x:
            double x_0 = 0;
            double x_np1 = 1;

            // Step value:
            double h = (x_np1-x_0)/(n+1);

            // Create xList, the list of all x values (n+2 points including end points):
            double *xList = new double[n+2];
            xList[0] = x_0;
            xList[n+1] = x_np1;
            for (int i=1; i<=n; i++) { xList[i] = x_0 + i*h; } 

            // Make the diagonal arrays for the general algorithm:
            double *a = new double[n-1];    // Lower diagonal.
            double *b = new double[n];      // Middle diagonal.
            double *c = new double[n-1];    // Upper diagonal.
            
            
            // General algorithm:           
            for (int i=0; i<=n-2; i++) { // Reset the diagonal elements.
                a[i] = -1;
                b[i] = 2;
                c[i] = -1;
            }
            b[n-1] = 2; // Fill the last element of b.

            // Get the time of the general algorithm:
            double tGen = general_algo(solutionGeneral, n, x_0, x_np1, a, b, c);
            tSumGen = tSumGen + tGen;  // Add to the sum (for calculating the average time).

            // Special algorithm:
            // Get the time of the special algorithm:
            double tSpec = special_algo(solutionSpecial, n, x_0, x_np1, a[0], b[0]);
            tSumSpec = tSumSpec + tSpec;

            cout << "n: " << n << ". tGen: " << tGen << ". tSpec: " << tSpec << "." << endl;

            // Delete allocated memory:
            delete[] a; delete[] b; delete[] c;
            delete[] xList; delete[] solutionGeneral; delete[] solutionSpecial; 
        
        }

        // Average times:
        double tAverageGen = tSumGen / ((double) numAverage);
        double tAverageSpec = tSumSpec / ((double) numAverage);
        timeListGeneral[p-1] = tAverageGen;
        timeListSpecial[p-1] = tAverageSpec;

        
    }

    // Write the time data to a .csv file.
    string fileNameGeneral = "timeGeneral.txt";
    string fileNameSpecial = "timeSpecial.txt";

    writeDataToCSV(nList, timeListGeneral, listLength, fileNameGeneral);
    writeDataToCSV(nList, timeListSpecial, listLength, fileNameSpecial);
    
    delete[] nList; delete[] timeListGeneral; delete[] timeListSpecial;
    return 0;
}
