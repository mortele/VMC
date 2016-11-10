#include <iostream>
#include <armadillo>
#include "system.h"
#include "Cores/atom.h"
#include "unittest.h"
#include "WaveFunctions/hydrogenwavefunction.h"
#include "WaveFunctions/directevaluationslater.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

int main(int, char**) {
    UnitTest::runAllTests();
}

