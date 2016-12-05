#include "unittest.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/gaussianslater.h"
#include <iostream>
#include <iomanip>
#include <mathimf.h>
#include <cmath>
#include <chrono>
#include "Math/exponentialapproximations.h"

using std::cout;
using std::endl;
using std::exp;
using std::log2;
using std::setprecision;



int main(int, char**) {
    UnitTest::runAllTests();




    /*
    double      dx   = 0.1;
    double      xMin = -500;
    double      xMax =  20;
    long int    N    = (xMax - xMin) / dx;

    double values  [11];
    double maxError[11];
    double xvalues [11];
    for (int i=0; i<11; i++) values  [i] = 0;
    for (int i=0; i<11; i++) maxError[i] = 0;

    for (long int i = 0; i < N; i++) {
        const double x   = xMin + N * dx;
        const double EXP = exp(x);

        values[2]  = expFast2 (x);
        values[3]  = expFast3 (x);
        values[4]  = expFast4 (x);
        values[5]  = expFast5 (x);
        values[6]  = expFast6 (x);
        values[7]  = expFast7 (x);
        values[8]  = expFast8 (x);
        values[9]  = expFast9 (x);
        values[10] = expFast10(x);

        for (int i = 2; i < 11; i++) {
            const double diff = fabs(EXP-values[i]);
            if (diff > maxError[i]) {
                maxError[i] = diff;
                xvalues [i] = x;
            }
        }
    }

    cout << "MAX ERRORS:" << endl;
    for (int i = 2; i < 11; i++) {
        cout << i << " : " << setprecision(20) << maxError[i] << "    x: " << xvalues[i] << endl;
    }


*/

    return 0;




/*
    const int N = 1e7;
    double sum = 0;
    double totalTime = 0;
    double bestTime  = 10;

    for (int k=0; k < 100; k++) {
        auto start = std::chrono::steady_clock::now();
        for (double i = 0.12534123; i < N; i++) {
            const double a = std::exp(-1.24124353242313 * i);
            sum += a;
        }
        auto finish = std::chrono::steady_clock::now();
        double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
        totalTime += elapsed_seconds;
        if (elapsed_seconds < bestTime) bestTime = elapsed_seconds;
    }
    std::cout << sum << std::endl;
    std::cout << std::setprecision(15) << totalTime << std::endl;
    std::cout << std::setprecision(15) << bestTime << std::endl;
    */
}

