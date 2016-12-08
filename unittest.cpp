#include "unittest.h"
#include "system.h"
#include "sampler.h"
#include "Cores/atom.h"
#include "Cores/harmonicoscillator.h"
#include "WaveFunctions/heliumwavefunction.h"
#include "WaveFunctions/hydrogenwavefunction.h"
#include "WaveFunctions/heliumwithjastrow.h"
#include "WaveFunctions/directevaluationslater.h"
#include "WaveFunctions/directevaluationslaterwithjastrow.h"
#include "WaveFunctions/gaussianslater.h"
#include "WaveFunctions/harmonicoscillatorwavefunction.h"
#include "RandomNumberGenerator/random.h"
#include "hartreefockbasisparser.h"
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
    //cout << "Running test: "; if (! testHydrogen())                          return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testNonInteractingHelium())              return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testHelium())                            return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testNumericalLaplacian())                return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testHeliumWithJastrowNumerical())        return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testDirectSlaterHelium())                return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testDirectSlaterWithJastrowHelium())     return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testDirectSlaterBeryllium())             return false; else cout << " -- passed" << endl;
    //cout << "Running test: "; if (! testDirectSlaterWithJastrowBeryllium())  return false; else cout << " -- passed" << endl;
    //---->cout << "Running test: "; if (! testGaussianSlaterHydrogenMolecule())    return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! HO3d())    return false; else cout << " -- passed" << endl;
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
                                                     nSpinDown));
    test->runMetropolisSilent(100000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = -4;
    assert(fabs(ref - E) < 1e-3);
    return true;
}

bool UnitTest::testDirectSlaterBeryllium() {
    printf("%-40s", "Direct eval. slater (Berillyium)"); fflush(stdout);
    int     nSpinUp     = 2;
    int     nSpinDown   = 2;
    double  alpha       = 4.0;
    System* test = setupNewTestSystem();
    test->setStepLength(1.0);
    test->setElectronInteraction(false);
    test->addCore(new Atom(test, zeros<vec>(3), 4));
    test->setWaveFunction(new DirectEvaluationSlater(test,
                                                     alpha,
                                                     nSpinUp,
                                                     nSpinDown));
    test->runMetropolis(100000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = -20;
    assert(fabs(ref - E) < 1e-3);
    return true;
}

bool UnitTest::testDirectSlaterWithJastrowHelium() {
    printf("%-40s", "Direct eval. slater /w Jastrow (Helium)"); fflush(stdout);
    int     nSpinUp     = 1;
    int     nSpinDown   = 1;
    double  alpha       = 1.843;
    double  beta        = 0.347;
    System* test = setupNewTestSystem();
    test->setElectronInteraction(true);
    test->addCore(new Atom(test, zeros<vec>(3), 2));
    test->setWaveFunction(new DirectEvaluationSlaterWithJastrow(test,
                                                                alpha,
                                                                beta,
                                                                nSpinUp,
                                                                nSpinDown));
    test->runMetropolisSilent(1000000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = -2.8901; // FYS4411 project 2
    assert(fabs(ref - E) < 1e-2);
    return true;
}

bool UnitTest::testDirectSlaterWithJastrowBeryllium()   {
    printf("%-40s", "Direct eval. slater /w Jastrow (Beryllium)"); fflush(stdout);
    int     nSpinUp     = 2;
    int     nSpinDown   = 2;
    double  alpha       = 3.983;
    double  beta        = 0.094;
    System* test = setupNewTestSystem();
    test->setElectronInteraction(true);
    test->setStepLength(1.5);
    test->addCore(new Atom(test, zeros<vec>(3), 4));
    test->setWaveFunction(new DirectEvaluationSlaterWithJastrow(test,
                                                                alpha,
                                                                beta,
                                                                nSpinUp,
                                                                nSpinDown));
    test->runMetropolis(1000000);
    const double E   = test->getSampler()->getEnergy();
    const double ref = -14.50; // FYS4411 project 2
    assert(fabs(ref - E) < 1e-2);
    return true;
}

bool UnitTest::testGaussianSlaterHydrogenMolecule() {
    printf("%-40s", "Slater (HF basis) (H molecule)"); fflush(stdout);
    System* test = setupNewTestSystem();
    HartreeFockBasisParser* parser = new HartreeFockBasisParser();
    parser->parseBasisFile("../../HartreeFock/HartreeFockBases/basis-2016-11-16-17.37.51");
    test->setWaveFunction(new GaussianSlater(test, parser));
    test->setStepLength(0.1);
    test->runMetropolis(1000000);
    return true;
}

bool UnitTest::HO3d() {
    printf("%-40s", "3d HO"); fflush(stdout);
    System* test = setupNewTestSystem();
    arma::vec pos = {0,0,0};
    double alpha = 1.1;
    double beta  = 0.56;
    double omega = 1.0;
    test->setElectronInteraction(true);
    test->setStepLength(1.5);
    test->setWaveFunction(new HarmonicOscillatorWaveFunction(test,
                                                             alpha,
                                                             beta,
                                                             omega));
    test->addCore(new HarmonicOscillator(test, pos, 2, omega));
    test->runMetropolis(5000000);
    return true;
}










