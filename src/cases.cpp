#include "cases.h"
#include "system.h"
#include "sampler.h"
#include "metropolis.h"
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




bool Cases::firstExample() {
    int     Z       = 2;
    double  alpha   = 1.843;
    double  beta    = 0.347;
    vec     position{0,0,0};

    System  He;
    //He.setStepLength(3.0);
    He.setImportanceSampling(true);
    He.addCore          (new Atom               (&He, position, Z));
    He.setWaveFunction  (new SlaterWithJastrow  (&He, beta,true));
    He.setOrbital       (new SlaterTypeOrbital  (alpha));
    He.runMetropolis    ((int) 1e7);
}

bool Cases::secondExample() {
    vec     position{0,0,0};
    System  Ne;
    Ne.setImportanceSampling    (true);
    Ne.setElectronInteraction   (false);
    Ne.setStepLength    (0.001);
    Ne.addCore          (new Atom               (&Ne, position, 10, 5, 4));
    Ne.setWaveFunction  (new SlaterWithJastrow  (&Ne, -1, false));
    Ne.setOrbital       (new HydrogenOrbital    (10.0));
    Ne.runMetropolis    ((int) 1e7);
}

bool Cases::thirdExample() {
    double  beta = 0.094;
    string  basisFileName = "Be-STO-3G";
    System  Be;
    Be.setImportanceSampling(true);
    Be.setWaveFunction  (new SlaterWithJastrow  (&Be, beta, true));
    Be.setOrbital       (new GaussianOrbital    (&Be, basisFileName));
    Be.runMetropolis    ((int) 1e7);
}

bool Cases::optimizeExample() {
    int     Z       = 2;
    double  alpha   = 1.843;
    double  beta    = 0.2;//0.347;
    vec     position{0,0,0};

    System  He;
    //He.setStepLength(3.0);
    He.setImportanceSampling(true);
    He.addCore          (new Atom               (&He, position, Z));
    He.setWaveFunction  (new SlaterWithJastrow  (&He, beta,true));
    He.setOrbital       (new SlaterTypeOrbital  (alpha));
    He.optimizeBeta(beta,0.001,50,1e6);
    He.getMetropolis()->runSteps((int) 1e7);
}

bool Cases::optimizeGaussianExample() {
    double beta = 0.36;
    System  He;
    He.setStepLength(0.08);
    He.setImportanceSampling(true);
    He.setWaveFunction  (new SlaterWithJastrow  (&He, beta, true));
    He.setOrbital       (new GaussianOrbital    (&He, "He-STO-6G-HF"));
    He.optimizeBeta(beta,0.001,50,1e6);
}

bool Cases::closedShellAtom() {
    boost::timer t;
    System* test = new System();
    test->setElectronInteraction(true);
    test->setImportanceSampling (true);
    test->setStepLength(0.025);
    //=========================================================================
    /*====================*/ string atom     = "He";
    /*====================*/ string orbital  = "Slater";
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

