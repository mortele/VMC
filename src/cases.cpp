#include "cases.h"
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


bool Cases::closedShellAtom() {
    boost::timer t;
    System* test = new System();
    test->setElectronInteraction(true);
    test->setImportanceSampling (true);
    test->setStepLength(0.025);
    //=========================================================================
    /*====================*/ string atom     = "Be";
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

