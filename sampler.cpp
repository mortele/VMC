#include "sampler.h"
#include "system.h"
#include "hamiltonian.h"
#include "electron.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setw;

void Sampler::sample(bool acceptedStep) {
    if (acceptedStep) {
        m_currentEnergy = m_hamiltonian->computeLocalEnergy();
    } else {
        m_currentEnergy = m_hamiltonian->getLocalEnergy();
    }
    m_currentEnergySquared            = m_currentEnergy * m_currentEnergy;
    m_cumulativeEnergy               += m_currentEnergy;
    m_cumulativeEnergySquared        += m_currentEnergySquared;
    m_blockCumulativeEnergy          += m_currentEnergy;
    m_blockCumulativeEnergySquared   += m_currentEnergySquared;
    m_blockAccepted                  += acceptedStep;
    m_totalAccepted                  += acceptedStep;
    m_totalSamplesTaken              += 1;
    m_blockSamplesTaken              += 1;
}

void Sampler::computeAverages() {
    const double E             = m_cumulativeEnergy        / m_totalSamplesTaken;
    const double E2            = m_cumulativeEnergySquared / m_totalSamplesTaken;
    m_acceptanceRate           = m_totalAccepted           / m_totalSamplesTaken;
    m_energy                   = E;
    m_variance                 = E2-E*E;
    m_numberOfMetropolisSteps  = m_totalSamplesTaken;
}

void Sampler::computeBlockAverages() {
    const double N        = m_blockSamplesTaken;
    const double E        = m_blockCumulativeEnergy        / N;
    const double E2       = m_blockCumulativeEnergySquared / N;
    m_blockEnergy         = E;
    m_blockVariance       = E2 - E*E;
    m_blockAcceptanceRate = m_blockAccepted / N;

    m_blockSamplesTaken            = 0;
    m_blockAccepted                = 0;
    m_blockCumulativeEnergy        = 0;
    m_blockCumulativeEnergySquared = 0;
}

Sampler::Sampler(System* system) {
    m_system = system;
}

void Sampler::setup() {
    m_hamiltonian = m_system->getHamiltonian();
}
