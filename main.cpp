#include "unittest.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/gaussianslater.h"
#include <iostream>
#include <iomanip>
//#include <mathimf.h>
#include <cmath>
#include <chrono>
#include "Math/exponentialapproximations.h"
#include "RandomNumberGenerator/random.h"

#ifdef ARMA_NO_DEBUG
    #undef ARMA_NO_DEBUG
#endif

// #define ARMA_NO_DEBUG

// Location of the armadillo 'config.hpp' file:
// /usr/local/include/armadillo_bits/config.hpp

using std::cout;
using std::endl;
using std::exp;
using std::log2;
using std::setprecision;

int main(int, char**) {
    UnitTest::runAllTests();
    return 0;
}

