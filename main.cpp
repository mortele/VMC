#include <iostream>
#include <armadillo>
#include "system.h"
#include "Cores/atom.h"
#include "WaveFunctions/hydrogenwavefunction.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

int main(int, char**) {
    System system;
    system.addCore(new Atom(&system, zeros<vec>(3), 1));
    system.setWaveFunction(new HydrogenWaveFunction(&system, 0.999));
    system.runMetropolis(3e6);
    return 0;
}
