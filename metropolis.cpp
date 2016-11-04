#include "metropolis.h"
#include "system.h"
#include "sampler.h"
#include "electron.h"
#include "Cores/core.h"
#include "WaveFunctions/wavefunction.h"
#include "RandomNumberGenerator/random.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::log10;
using std::pow;

bool Metropolis::step() {
    int    electron       = Random::nextInt(0, m_numberOfElectrons  - 1);
    int    dimension      = Random::nextInt(0, m_numberOfDimensions - 1);
    double proposedChange = Random::nextDouble( - m_stepLengthHalf,
                                                  m_stepLengthHalf);
    m_system->getElectrons().at(electron)->adjustPosition(proposedChange,
                                                          dimension);
    m_waveFunction->updateOldWaveFunctionValue();
    const double R = m_waveFunction->computeWaveFunctionRatio(electron);

    if (R > Random::nextDouble(0, 1)) {
        return true;
    } else {
        m_system->getElectrons().at(electron)->adjustPosition(- proposedChange,
                                                                dimension);
        return false;
    }
}

Metropolis::Metropolis(System* system) {
    m_system = system;
    m_stepLengthHalf = 0.5 * m_stepLength;
}

void Metropolis::setup() {
    m_sampler            = m_system->getSampler();
    m_waveFunction       = m_system->getWaveFunction();
    m_numberOfElectrons  = m_system->getNumberOfElectrons();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
}

void Metropolis::runSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
    printInitialInfo();
    m_waveFunction->evaluateWaveFunctionInitial();
    for (int i = 0; i < steps; i++) {
        bool acceptedStep = step();
        m_sampler->sample(acceptedStep);
        printIterationInfo(i);
    }
    m_sampler->computeAverages(steps);
    printFinalInfo();
}

void Metropolis::printInitialInfo() {
    printf(" =============== Starting Metropolis Algorithm ============== \n");
    printf(" => Number of steps:       %-10d \n", (int) m_numberOfMetropolisSteps);
    printf(" => Number of dimensions:  %-10d \n", (int) m_numberOfDimensions);
    printf(" => Number of electrons:   %-10d \n", (int) m_numberOfElectrons);
    printf(" => Number of cores:       %-10d \n", (int) m_system->getCores().size());
    printf("      ------------------------------------------------------- \n");
    for (Core* core : m_system->getCores()) {
        printf("      | %-20s (%5.3f, %5.3f, %5.3f)          | \n", core->getInfo().c_str(),
                                                                    core->getPosition()(0),
                                                                    core->getPosition()(1),
                                                                    core->getPosition()(2));
    }
    printf("      ------------------------------------------------------- \n\n");
    printf(" ============================================================ \n");
    printf(" %17s %5s %14s %10s \n", " ", "Total", " ", "Block" );
    printf(" %5s %12s %12s %12s %12s\n", "Step", "Energy", "Variance", "Energy", "Variance");
    printf(" ------------------------------------------------------------ \n");
    fflush(stdout);
}

void Metropolis::printIterationInfo(int iteration) {
    const int skip = 100000;
    if (iteration != 0 && iteration % (20 * skip) == 0) {
        printf(" ============================================================ \n");
        printf(" %17s %5s %14s %10s \n", " ", "Total", " ", "Block" );
        printf(" %5s %12s %12s %12s %12s\n", "Step", "Energy", "Variance", "Energy", "Variance");
        printf(" ------------------------------------------------------------ \n");
        fflush(stdout);
    }
    if (iteration != 0 && iteration % skip == 0) {
        const int exponent  = log10(iteration);
        const int preFactor = iteration / pow(10, exponent);

        printf(" %3de%2-d %12.5g %12.5g %12.5g %12.5g \n",
                preFactor,
                exponent,
                m_sampler->getEnergyAverage(iteration),
                m_sampler->getVariance(iteration),
                m_sampler->getEnergyBlockAverage(iteration),
                m_sampler->getVarianceBlock(iteration));
    }

}

void Metropolis::printFinalInfo() {
    printf(" ============================================================ \n");
    printf("\n Self consistency SUCCESFULLY reached. \n\n");
    printf(" => Metropolis steps:       %30d   \n",  m_numberOfMetropolisSteps);
    printf(" => Final energy average:   %30.16g \n", m_sampler->getEnergy());
    printf(" => Final variance:         %30.16g \n", m_sampler->getVariance());
    printf(" ============================================================ \n");
    fflush(stdout);
}
















