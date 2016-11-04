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
    system.setWaveFunction(new HydrogenWaveFunction(&system));
    system.runMetropolis(5);

    cout << "Hello World!" << endl;
    return 0;
}
