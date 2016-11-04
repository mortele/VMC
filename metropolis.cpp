#include "metropolis.h"
#include "system.h"
#include "sampler.h"
#include "electron.h"
#include "WaveFunctions/wavefunction.h"
#include "RandomNumberGenerator/random.h"
#include <iostream>

using std::cout;
using std::endl;

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
    m_waveFunction->evaluateWaveFunctionInitial();
    for (int i = 0; i < steps; i++) {
        bool acceptedStep = step();
        m_sampler->sample(acceptedStep);
    }
}
