#include "unittest.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/gaussianslater.h"
#include <iostream>
#include <iomanip>
//#include <mathimf.h>
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
    return 0;
}

