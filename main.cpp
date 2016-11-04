#include <iostream>
#include <armadillo>
#include "system.h"
#include "Cores/atom.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

int main(int, char**) {
    System system;
    Atom atom(&system, zeros<vec>(3), 1);
    system.addCore(&atom);

    cout << "Hello World!" << endl;
    return 0;
}
