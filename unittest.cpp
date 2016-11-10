#include "unittest.h"
#include "system.h"
#include "sampler.h"
#include "Cores/atom.h"
#include "WaveFunctions/heliumwavefunction.h"
#include "WaveFunctions/hydrogenwavefunction.h"
#include "WaveFunctions/heliumwithjastrow.h"
#include "WaveFunctions/directevaluationslater.h"
#include "RandomNumberGenerator/random.h"
#include <armadillo>
#include <cassert>
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;
using arma::zeros;
using arma::vec;

System* UnitTest::setupNewTestSystem() {
    Random::seed(0);
    System* test = new System();
    return test;
}

bool UnitTest::runAllTests() {
    cout << "Running all tests." << endl;
    cout << "=================================================================" << endl;
    //cout << "Running test: "; if (! testHydrogen())                     return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testNonInteractingHelium())         return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testHelium())                       return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testNumericalLaplacian())           return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testHeliumWithJastrowNumerical())   return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testDirectSlaterHelium())           return false; else cout << " -- passed" << endl;
    cout << "=================================================================" << endl;
    cout << "All tests passed." << endl;
    return true;
}

bool UnitTest::testHydrogen() {
    printf("%-40s", "Hydrogen"); fflush(stdout);
    System* test = setupNewTestSystem();
    test->setStepLength(5.0);
    test->addCore(new Atom(test, zeros<vec>(3), 1));
    test->setWaveFunction(new HydrogenWaveFunction(test, 1.0));
    test->runMetropolisSilent(100000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = -0.5;
    assert(E == ref);
    return true;
}

bool UnitTest::testNonInteractingHelium() {
    printf("%-40s", "Non-interacting Helium"); fflush(stdout);
    System* test = setupNewTestSystem();
    test->setStepLength(3.0);
    test->setElectronInteraction(false);
    test->addCore(new Atom(test, zeros<vec>(3), 2));
    test->setWaveFunction(new HeliumWaveFunction(test, 2.0));
    test->runMetropolisSilent(100000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = -4.0;
    assert(E == ref);
    return true;

}

bool UnitTest::testNumericalLaplacian() {
    printf("%-40s", "Numerical laplacian (Non-Int. Helium)"); fflush(stdout);
    System* test = setupNewTestSystem();
    //test->setStepLength(1.17163);
    test->setElectronInteraction(false);
    test->addCore(new Atom(test, zeros<vec>(3), 2));
    test->setWaveFunction(new HeliumWaveFunction(test, 2.0, true));
    test->getWaveFunction()->setStepLength(1e-5);
    test->runMetropolisSilent(100000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = -4.0;
    assert(fabs(ref - E) < 1e-3);
    return true;
}

bool UnitTest::testDirectSlaterHelium() {
    printf("%-40s", "Direct eval. slater (Helium)"); fflush(stdout);
    int     nSpinUp     = 1;
    int     nSpinDown   = 1;
    double  alpha       = 2.0;
    System* test = setupNewTestSystem();
    test->setElectronInteraction(false);
    test->addCore(new Atom(test, zeros<vec>(3), 2));
    test->setWaveFunction(new DirectEvaluationSlater(test,
                                                     alpha,
                                                     nSpinUp,
                                                     nSpinDown,
                                                     true));
    test->runMetropolis(100000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = -4;
    cout << E   << endl;
    cout << ref << endl;
    assert(fabs(ref - E) < 1e-3);
    return true;
}

bool UnitTest::testHelium() {
    printf("%-40s", "Helium"); fflush(stdout);
    System* test = setupNewTestSystem();
    test->addCore(new Atom(test, zeros<vec>(3), 2));
    test->setWaveFunction(new HeliumWaveFunction(test, 27.0 / 16.0));
    test->setStepLength(12.0);
    test->runMetropolisSilent(1000000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = 0.5*std::pow(3.0/2.0,6)*(-0.5); // Griffiths pp. 303
    assert(fabs(ref - E) < 3e-3);
    return true;
}

bool UnitTest::testHeliumWithJastrowNumerical() {
    printf("%-40s", "Helium with Jastrow (numerical)"); fflush(stdout);
    System* test = setupNewTestSystem();
    test->addCore(new Atom(test, zeros<vec>(3), 2));
    test->setWaveFunction(new HeliumWithJastrow(test, 1.843, 0.347, true));
    test->setStepLength(1.0);
    test->runMetropolisSilent(100000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = -2.8901; // FYS4411 project 2
    assert(fabs(ref - E) < 1e-2);
    return true;
}







