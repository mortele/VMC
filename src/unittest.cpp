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
#include "WaveFunctions/slaterwithjastrow.h"
#include "WaveFunctions/Orbitals/hydrogenorbital.h"
#include "WaveFunctions/Orbitals/slatertypeorbital.h"
#include "WaveFunctions/Orbitals/gaussianorbital.h"
#include "RandomNumberGenerator/random.h"
#include "hartreefockbasisparser.h"
#include <armadillo>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <boost/timer.hpp>
#include <fstream>
#include <string>


using std::cout;
using std::endl;
using std::setprecision;
using std::string;
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
    cout << "Running test: "; if (! testHydrogen())                          return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testNonInteractingHelium())              return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testHelium())                            return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testNumericalLaplacian())                return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testHeliumWithJastrowNumerical())        return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testDirectSlaterHelium())                return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testDirectSlaterWithJastrowHelium())     return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testDirectSlaterBeryllium())             return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testDirectSlaterWithJastrowBeryllium())  return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testGaussianSlaterHydrogenMolecule())    return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! HO3d())                                  return false; else cout << " -- passed" << endl;
    cout << "Running test: "; if (! testImportanceSampledSlaterWithJastrowBe())                                  return false; else cout << " -- passed" << endl;
    cout << testSlaterWithJastrowGaussian() << endl;
    cout << testSlaterWithJastrowGaussianHe() << endl;
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
    test->runMetropolis(100000);
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
    //parser->parseBasisFile("../HartreeFock/HartreeFockBases/basis-2017-08-31-20.29.20");
    parser->parseBasisFile(test, "../HartreeFock/HartreeFockBases/basis-2017-09-01-12.34.04");
    test->setWaveFunction(new GaussianSlater(test, parser));
    //test->setStepLength(2.5);
    test->runMetropolis((int) 1e7);
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

bool UnitTest::testImportanceSampledSlaterWithJastrowBe() {
    boost::timer t;

    Random::randomSeed();
    printf("%-40s", "Imp. sampling Slater w. Jastrow (Be)"); fflush(stdout);
    System* test = setupNewTestSystem();
    arma::vec pos = {0,0,0};
    double alpha = 3.983; // 1.843;
    double beta  = 0.094; // 0.347;
    test->setElectronInteraction(true);
    test->setImportanceSampling (true);
    test->setStepLength(0.01);
    test->setWaveFunction(new SlaterWithJastrow(test,beta,true));
    test->setOrbital     (new HydrogenOrbital(alpha));
    test->addCore        (new Atom(test,pos,4,2,2));
    test->runMetropolis((int) 1e6);

    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
    return true;
}

bool UnitTest::testSlaterWithJastrowGaussian() {
    boost::timer t;
    Random::randomSeed();


    //printf("%-40s", "Imp. sampling Slater w. Jastrow (Gaussian basis) (Be)"); fflush(stdout);

    std::ofstream outFile;
    outFile.open("betaEnergy.txt", std::ios::out);
    if (! outFile.is_open()) {
        cout << "coudlnt open file." << endl;
        exit(1);
    }

    int    n     = 100;
    double beta0 = 0.1;
    double beta1 = 0.9;
    double db    = (beta1-beta0)/n;

    for (int i = 0; i < n; i++) {
        System* test = setupNewTestSystem();
        arma::vec pos = {0,0,0};
        //double alpha = 1.843; // 3.983;
        double beta  = beta0+db*i; //0.4; //0.347; // 0.094;
        test->setElectronInteraction(true);
        test->setImportanceSampling (true);
        test->setStepLength(0.01);
        test->setWaveFunction(new SlaterWithJastrow(test,beta,true));
        test->setOrbital     (new GaussianOrbital(test, "He-321G"));
        //test->setOrbital     (new GaussianOrbital("He-6311++G**"));
        test->addCore        (new Atom(test,pos,2,1,1));
        double E = test->runMetropolisSilent((int) 5e5);
        outFile << setprecision(16) << beta << " " << E << endl;
        printf("%10g %10g \n", beta, E);
        if (i % 1 == 0) fflush(stdout);
    }
    outFile.close();
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
    return true;
}

bool UnitTest::testSlaterWithJastrowGaussianHe() {
    boost::timer t;
    Random::randomSeed();
    //Random::seed(925733851);
    System* test = setupNewTestSystem();
    test->setElectronInteraction(true);
    test->setImportanceSampling (true);
    test->setStepLength(0.025);
    //=========================================================================
    /*====================*/ string atom     = "Be";
    /*====================*/ string orbital  = "Gaussian";
    /*====================*/ string basis    = "Be-STO-6G";
    //=========================================================================
    double alpha, beta;
    if (atom=="He") {
        alpha = 1.843;
        beta  = 0.347;
        if (orbital != "Gaussian") test->addCore(new Atom(test,vec{0,0,0},2,1,1));
    } else if (atom == "Be") {
        alpha = 3.983;
        beta  = 0.094;
        if (orbital != "Gaussian") test->addCore(new Atom(test,vec{0,0,0},4,2,2));
    } else {
        alpha = 10.22;
        beta  = 0.091;
        if (orbital != "Gaussian") test->addCore(new Atom(test,vec{0,0,0},10,5,5));
    }
    test->setWaveFunction(new SlaterWithJastrow(test,beta,true));
    if (orbital == "Slater") {
        test->setOrbital(new SlaterTypeOrbital(alpha));
    } else if (orbital == "Hydrogen") {
        test->setOrbital(new HydrogenOrbital(alpha));
    } else {
        test->setOrbital(new GaussianOrbital(test, basis));
    }
    test->runMetropolis((int) 5e6);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
    return true;
}


// clang++ plain :          109.777
// icpc fast :              111.038
// clang++ -ffast-math :    110.777
// icpc -fp-model fast=2 :  107.399




