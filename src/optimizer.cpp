#include "optimizer.h"
#include "system.h"
#include "metropolis.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;



void Optimizer::setup() {
    m_waveFunction  = m_system->getWaveFunction();
    m_sampler       = m_system->getSampler();
    m_metropolis    = m_system->getMetropolis();

    if (m_waveFunction->containsJastrow()) {
        m_jastrow = true;
    }
}

double Optimizer::optimizeBeta(double   beta,
                               double   tollerance,
                               int      maxIterations,
                               int      cycles) {

    m_maximumIterations = maxIterations;
    m_cycles            = cycles;
    m_tollerance        = tollerance;
    m_beta              = beta;

    double oldBeta       = 0;
    double oldBetaChange = 0;
    double oldGradient   = 0;
    double gradient      = 0;
    double step          = 1;

    printInitialInfo();
    for (int iteration = 0; iteration < m_maximumIterations; iteration++) {
        m_system->getWaveFunction()->setBeta(m_beta);

        if (iteration == 0) m_system->runMetropolisSilent(cycles);
        else m_metropolis->runStepsSilent(cycles);

        // Gradient descent.
        oldGradient = gradient;
        oldBeta     = m_beta;
        gradient    = m_sampler->getJastrowGradient();
        double gradientChange = gradient - oldGradient;

        // Barzilai-Borwein method.
        if (iteration == 0) {
            step = 1;
        } else {
            step = oldBetaChange / gradientChange;
            step = 1;
        }
        printIterationInfo(iteration);
        m_beta -= step * gradient;

        m_sampler->clear();

        if (iteration != 0) {
            if (fabs(oldBeta - m_beta) < m_tollerance) {
                break;
            }
        }
        oldBetaChange = m_beta - oldBeta;
        oldGradient = gradient;
    }
    return m_beta;
}

Optimizer::Optimizer(System* system) {
    m_system = system;
}


void Optimizer::printInitialInfo() {
    printf(" =============== Starting Jastrow Optimization ============== \n");
    printf(" => Max iterations:        %-10g \n", m_maximumIterations);
    printf(" => Starting beta:         %-10g \n", m_beta);
    printf(" => Cycles per iteration:  %-10g \n", m_cycles);
    printf(" => Tollerance:            %-10g \n", m_tollerance);
    cout << endl;
    printf(" ==================================================================================== \n");
    printf(" %5s %12s %12s %12s %12s %12s \n",
           "Step", "Energy", "Std.dev.", "Acc. rate", "Beta", "Beta grad.");
    printf(" ------------------------------------------------------------------------------------ \n");

    fflush(stdout);
}

void Optimizer::printIterationInfo(int iteration) {
    printf(" %5d %12.5g %12.5g %12.5g %12.5g %12.5g \n",
           iteration,
           m_sampler->getEnergy(),
           m_sampler->getStandardDeviation(),
           m_sampler->getAcceptanceRate(),
           m_beta,
           m_sampler->getJastrowGradient());
    fflush(stdout);
}
