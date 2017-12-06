#include <iostream>
#include <iomanip>
#include "cases.h"

#ifdef ARMA_NO_DEBUG
    #undef ARMA_NO_DEBUG
#endif

// #define ARMA_NO_DEBUG

// Location of the armadillo 'config.hpp' file:
// /usr/local/include/armadillo_bits/config.hpp


int main(int, char**) {
    //Cases::firstExample();
    //Cases::secondExample();
    //Cases::thirdExample();
    //Cases::optimizeExample();
    //Cases::optimizeGaussianExample();
    Cases::nonInteracting();
    //Cases::closedShellAtom();
    return 0;
}

