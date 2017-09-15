#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <armadillo>
#include "system.h"

#include "Cores/atom.h"
#include "Cores/harmonicoscillator.h"

#include "WaveFunctions/hydrogenwavefunction.h"
#include "WaveFunctions/heliumwavefunction.h"
#include "WaveFunctions/slaterwithjastrow.h"

#include "WaveFunctions/Orbitals/slatertypeorbital.h"
#include "WaveFunctions/Orbitals/gaussianorbital.h"
#include "WaveFunctions/Orbitals/hydrogenorbital.h"


using namespace arma;
using namespace std;

TEST_CASE( " s ", "[Slater]") {
    System* test = new System();
    test->setStepLength(1.0);
    test->addCore(new Atom(test, zeros<vec>(3), 1));
    test->setWaveFunction(new HydrogenWaveFunction(test, 1.0));

    REQUIRE( test->runMetropolisSilent(10000) == Approx(-0.5) );
}

TEST_CASE(" a ", "[Slater]") {
    System* test = new System();
    test->setStepLength(1.0);
    test->setElectronInteraction(false);
    test->addCore(new Atom(test, zeros<vec>(3), 2));
    test->setWaveFunction(new HeliumWaveFunction(test, 2.0));

    REQUIRE( test->runMetropolisSilent(10000) == Approx(-4.0) );
}
