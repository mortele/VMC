#include "unittest.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/gaussianslater.h"
#include <iostream>
#include <iomanip>
//#include <mathimf.h>
#include <cmath>
#include <chrono>
#include "Math/exponentialapproximations.h"

#ifdef ARMA_NO_DEBUG
    #undef ARMA_NO_DEBUG
#endif

// #define ARMA_NO_DEBUG

// Location of the armadillo 'config.hpp' file:
// /usr/local/include/armadillo_bits/config.hpp

using std::cout;
using std::endl;
using std::setprecision;


int main(int, char**) {
    cout << std::numeric_limits<float>::epsilon()*100 << endl;
    //STOFitter fit(3.98);
    //fit.computeFit(3);

    //UnitTest::runAllTests();
    return 0;
}

