#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <armadillo>
#include "system.h"
#include "sampler.h"

#include "Cores/atom.h"
#include "Cores/harmonicoscillator.h"

#include "WaveFunctions/hydrogenwavefunction.h"
#include "WaveFunctions/heliumwavefunction.h"
#include "WaveFunctions/slaterwithjastrow.h"
#include "WaveFunctions/harmonicoscillatorwavefunction.h"

#include "WaveFunctions/Orbitals/slatertypeorbital.h"
#include "WaveFunctions/Orbitals/gaussianorbital.h"
#include "WaveFunctions/Orbitals/hydrogenorbital.h"


using namespace arma;
using namespace std;

const int monteCarloCycles = (int) 1e3;

TEST_CASE("Simple hydrogen", "[Hydrogen wave functions]") {
    System* test = new System();
    test->setStepLength(2.0);
    test->addCore(new Atom(test, vec{0,0,0}, 1));
    test->setWaveFunction(new HydrogenWaveFunction(test, 1.0));

    REQUIRE( test->runMetropolisSilent(monteCarloCycles) == Approx(-0.5) );
    REQUIRE( test->getSampler()->getVariance() == Approx(0.0) );
}

TEST_CASE("Non-interacting Helium ", "[Hydrogen wave functions]") {
    System* test = new System();
    test->setStepLength(2.0);
    test->setElectronInteraction(false);
    test->addCore(new Atom(test, vec{0,0,0}, 2));
    test->setWaveFunction(new HeliumWaveFunction(test, 2.0));

    REQUIRE( test->runMetropolisSilent(monteCarloCycles) == Approx(-4.0) );
    REQUIRE( test->getSampler()->getVariance() == Approx(0.0) );
}

TEST_CASE("Non-interacting 3D quantum dot", "[Quantum dot]") {
    System* test = new System();
    test->setStepLength(2.0);
    test->setElectronInteraction(false);
    test->addCore(new HarmonicOscillator(test, vec{0,0,0}, 2, 1.0));
    test->setWaveFunction(new HarmonicOscillatorWaveFunction(test,1.0, 0, 1.0));

    REQUIRE( test->runMetropolisSilent(monteCarloCycles) == Approx(3.0) );
    REQUIRE( test->getSampler()->getVariance() == Approx(0.0) );
}

TEST_CASE("Non-interacting Helium", "[Slater]") {
    System* test = new System();
    test->setStepLength(2.0);
    test->setElectronInteraction(false);
    test->addCore(new Atom(test, vec{0,0,0}, 2));
    test->setWaveFunction(new SlaterWithJastrow(test,0,false));
    test->setOrbital(new HydrogenOrbital(2.0));

    REQUIRE( test->runMetropolisSilent(monteCarloCycles) == Approx(-4.0) );
    REQUIRE( test->getSampler()->getVariance() == Approx(0.0) );
}

TEST_CASE("Non-interacting Beryllium", "[Slater]") {
    System* test = new System();
    test->setStepLength(2.0);
    test->setElectronInteraction(false);
    test->addCore(new Atom(test, vec{0,0,0}, 4));
    test->setWaveFunction(new SlaterWithJastrow(test,0,false));
    test->setOrbital(new HydrogenOrbital(4.0));

    REQUIRE( test->runMetropolisSilent(monteCarloCycles) == Approx(-20.0) );
    REQUIRE( test->getSampler()->getVariance() == Approx(0.0) );
}

TEST_CASE("Non-interacting Neon", "[Slater]") {
    System* test = new System();
    test->setStepLength(2.0);
    test->setElectronInteraction(false);
    test->addCore(new Atom(test, vec{0,0,0}, 10));
    test->setWaveFunction(new SlaterWithJastrow(test,0,false));
    test->setOrbital(new HydrogenOrbital(10.0));

    REQUIRE( test->runMetropolisSilent(monteCarloCycles) == Approx(-200.0) );
    REQUIRE( test->getSampler()->getVariance() == Approx(0.0) );
}

TEST_CASE("Gaussian Be STO-6G orbitals", "[Gaussian orbitals]") {
    System* test2 = new System();
    test2->addCore(new Atom(test2, vec{0,0,0}, 4));
    test2->setWaveFunction(new SlaterWithJastrow(test2,0,false));

    GaussianOrbital   gaussian = GaussianOrbital(test2,"Be-STO-6G");
    SlaterTypeOrbital slater   = SlaterTypeOrbital(3.983);
    double x =  1.2345;
    double y = -0.2872;
    double z =  0.9258;
    double STO   = 1; //slater.evaluate(x,y,z,0,0);
    double STO6G = gaussian.evaluate(x,y,z,0,0);

    cout << STO << " " << STO6G << endl;
    REQUIRE( STO == Approx(STO6G) );
}


