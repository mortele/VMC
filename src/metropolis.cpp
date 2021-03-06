#include "metropolis.h"
#include "system.h"
#include "sampler.h"
#include "electron.h"
#include "hamiltonian.h"
#include "Cores/core.h"
#include "WaveFunctions/wavefunction.h"
#include <iostream>
#include <cmath>
#include <iomanip>

using std::cout;
using std::endl;
using std::log10;
using std::pow;

bool Metropolis::step() {
    std::uniform_int_distribution<int> uniformIntElectron {0, m_numberOfElectrons-1};
    int    electron       = uniformIntElectron (m_randomGenerator);
    int    dimension      = -1;
    double proposedChange;
    double xProposedChangeImportanceSampling;
    double yProposedChangeImportanceSampling;
    double zProposedChangeImportanceSampling;
    if (! m_importanceSampling) {
        std::uniform_int_distribution<int> uniformIntDimension{0, m_numberOfDimensions-1};
        dimension      = uniformIntDimension(m_randomGenerator);
        std::uniform_real_distribution<double> uniformDouble{-m_stepLengthHalf,m_stepLengthHalf};
        proposedChange = uniformDouble(m_randomGenerator);
    } else {
        double D = 0.5;

        std::normal_distribution<double> normalDistribution{0,1};

        xProposedChangeImportanceSampling
                = normalDistribution(m_randomGenerator) *
                  m_dtSqrt + m_waveFunction->getQuantumForceOld(electron,0) *
                  m_dt * D;

        yProposedChangeImportanceSampling
                = normalDistribution(m_randomGenerator) *
                  m_dtSqrt + m_waveFunction->getQuantumForceOld(electron,1) *
                  m_dt * D;

        zProposedChangeImportanceSampling
                = normalDistribution(m_randomGenerator) *
                  m_dtSqrt + m_waveFunction->getQuantumForceOld(electron,2) *
                  m_dt * D;
    }

    m_waveFunction->passProposedChangeToWaveFunction(electron, dimension);
    if (! m_importanceSampling) {
        m_system->getElectrons().at(electron)->adjustPosition(proposedChange, dimension);
    } else {
        m_system->getElectrons().at(electron)->adjustPosition(xProposedChangeImportanceSampling, 0);
        m_system->getElectrons().at(electron)->adjustPosition(yProposedChangeImportanceSampling, 1);
        m_system->getElectrons().at(electron)->adjustPosition(zProposedChangeImportanceSampling, 2);
    }
    m_waveFunction->updateOldWaveFunctionValue();
    double R = m_waveFunction->computeWaveFunctionRatio(electron);
    if (m_importanceSampling) {
        double greensFunction = computeGreensFunction();
        R *= R * greensFunction;
    }

    std::uniform_real_distribution<double> uniformDouble{0,1};
    if (R > uniformDouble(m_randomGenerator)) {
        m_waveFunction->updateWaveFunctionAfterAcceptedStep();
        return true;
    } else {
        m_waveFunction->updateWaveFunctionAfterRejectedStep();
        if (! m_importanceSampling) {
            m_system->getElectrons().at(electron)->adjustPosition(-proposedChange, dimension);
        } else {
            m_system->getElectrons().at(electron)->adjustPosition(-xProposedChangeImportanceSampling, 0);
            m_system->getElectrons().at(electron)->adjustPosition(-yProposedChangeImportanceSampling, 1);
            m_system->getElectrons().at(electron)->adjustPosition(-zProposedChangeImportanceSampling, 2);
        }
        return false;
    }
}


Metropolis::Metropolis(System* system) {
    m_system = system;
    m_stepLengthHalf = 0.5 * m_stepLength;
}

void Metropolis::setImportanceSampling(bool importanceSampling) {
    m_importanceSampling = importanceSampling;
}

void Metropolis::setup() {
    m_sampler            = m_system->getSampler();
    m_waveFunction       = m_system->getWaveFunction();
    m_numberOfElectrons  = m_system->getNumberOfElectrons();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    if (! m_stepLengthSetManually) {
        if (! m_importanceSampling) {
            double minimumSize = 10.0;
            for (Core* core : m_system->getCores()) {
                minimumSize = std::min(minimumSize, core->getSize());
            }
            m_stepLength     = 4.0 * minimumSize;
            m_stepLengthHalf = 0.5 * m_stepLength;
        } else {
            m_dt        = 0.01;
            m_dtSqrt    = sqrt(m_dt);
        }
    }
}

void Metropolis::setStepLength(double stepLength) {
    m_stepLength     = stepLength;
    m_stepLengthHalf = 0.5 * m_stepLength;
    m_dt             = stepLength;
    m_dtSqrt         = sqrt(m_dt);
    m_stepLengthSetManually = true;
}

double Metropolis::runSteps(int steps) {
    setup();
    m_system->getHamiltonian()->computeCoreCorePotentialEnergy();
    m_system->getHamiltonian()->computeElectronCorePotentialEnergy();
    m_system->getHamiltonian()->computeElectronElectronPotentialEnergy();
    m_numberOfMetropolisSteps = steps;
    if (! m_silent) printInitialInfo();
    m_waveFunction->evaluateWaveFunctionInitial();
    for (m_i = 0; m_i < steps; m_i++) {
        bool acceptedStep = step();
        if (m_thermalize && (m_i > 2500*4)) m_sampler->sample(acceptedStep);
        printIterationInfo(m_i);
    }
    m_sampler->computeAverages();
    if (! m_silent) printFinalInfo();
    return m_sampler->getEnergy();
}

double Metropolis::runStepsSilent(int steps) {
    m_silent = true;
    double E = runSteps(steps);
    m_silent = false;
    return E;
}

void Metropolis::printInitialInfo() {
    printf(" =============== Starting Metropolis Algorithm ============== \n");
    printf(" => Number of steps:       %-10g \n", (double) m_numberOfMetropolisSteps);
    printf(" => Number of dimensions:  %-10d \n", (int) m_numberOfDimensions);
    printf(" => Number of electrons:   %-10d \n", (int) m_numberOfElectrons);
    printf(" => Step length:           %-10g \n", m_stepLength);
    printf(" => Number of cores:       %-10d \n", (int) m_system->getCores().size());
    printf("      ------------------------------------------------------- \n");
    for (Core* core : m_system->getCores()) {
        printf("      | %-20s (%5.3f, %5.3f, %5.3f)          | \n", core->getInfo().c_str(),
               core->getPosition()(0),
               core->getPosition()(1),
               core->getPosition()(2));
    }
    printf("      ------------------------------------------------------- \n\n");
    printf(" ==================================================================================== \n");
    printf(" %18s %5s %27s %10s \n", " ", "Total", " ", "Block" );
    printf(" %5s %12s %12s %12s %12s %12s %12s\n", "Step", "Energy", "Std.dev.", "Energy", "Variance", "Acc. rate", "Virial r.");
    printf(" ------------------------------------------------------------------------------------ \n");
    fflush(stdout);
}

void Metropolis::printIterationInfo(int iteration) {
    const int skip = 2500;
    if (! m_silent) {
        if (iteration != 0 && iteration % (20*20 * skip) == 0) {
            printf(" ==================================================================================== \n");
            printf(" %18s %5s %27s %10s \n", " ", "Total", " ", "Block" );
            printf(" %5s %12s %12s %12s %12s %12s %12s\n", "Step", "Energy", "Std.dev.", "Energy", "Variance", "Acc. rate", "Virial r.");
            printf(" ------------------------------------------------------------------------------------ \n");
        }
    }
    if (iteration != 0 && iteration % skip == 0) {
        const int exponent  = log10(iteration);
        const int preFactor = iteration / pow(10, exponent);

        m_sampler->computeBlockAverages();
        m_sampler->computeAverages();
        if (! m_silent) {
            if (iteration % (20 * skip) == 0) {
            printf(" %3de%-2d %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",
                   preFactor,
                   exponent,
                   m_sampler->getEnergy(),
                   m_sampler->getStandardDeviation(),
                   m_sampler->getBlockEnergy(),
                   m_sampler->getBlockVariance(),
                   m_sampler->getBlockAcceptanceRate(),
                   m_sampler->getBlockVirialRatio());
            }
        }
    }
    if (! m_silent) fflush(stdout);
}

void Metropolis::printFinalInfo() {
    if (! m_silent) {
        printf(" ======================================================================== \n");
        printf("\n Metropolis algorithm finished. \n\n");
        printf(" => Metropolis steps:       %30g   \n",  (double) m_numberOfMetropolisSteps);
        printf(" => Final acceptance rate:  %30.16g \n", m_sampler->getAcceptanceRate());
        printf(" => Final energy average:   %30.16g \n", m_sampler->getEnergy());
        printf(" => Final variance:         %30.16g \n", m_sampler->getVarianceBlocking());
        printf(" ============================================================ \n");
        fflush(stdout);
    }
}

double Metropolis::computeGreensFunction() {
    const double D = 0.5;
    WaveFunction* wf = m_waveFunction;
    double greensFunction = 0;
    for (int i = 0; i < m_numberOfElectrons; i++) {
        for (int j = 0; j < 3; j++) {
            greensFunction += 0.5 *
                 (wf->getQuantumForceOld(i,j) + wf->getQuantumForce(i,j)) *
                 (D * m_dt * 0.5 *
                    (wf->getQuantumForceOld(i,j) - wf->getQuantumForce(i,j)) -
                     wf->getPosition(i,j)        + wf->getPositionOld(i,j));
        }
    }
    return exp(greensFunction);
}












