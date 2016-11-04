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
    //cout << "[ "  << setw(10) << m_system->getElectrons().at(0)->getPosition()(0) << ", "
    //              << setw(10) << m_system->getElectrons().at(0)->getPosition()(1) << ", "
    //              << setw(10) << m_system->getElectrons().at(0)->getPosition()(2)
    //     << " ] " << setw(10) << m_currentEnergy << " : " << acceptedStep << endl;
    m_currentEnergySquared            = m_currentEnergy * m_currentEnergy;
    m_cumulativeEnergy               += m_currentEnergy;
    m_cumulativeEnergySquared        += m_currentEnergySquared;
    m_cumulativeBlockEnergy          += m_currentEnergy;
    m_cumulativeBlockEnergySquared   += m_currentEnergySquared;
    m_totalSamplesTaken              += 1;
}

void Sampler::computeAverages(int steps) {
    m_numberOfMetropolisSteps  = m_totalSamplesTaken;
    m_cumulativeEnergy        /= m_totalSamplesTaken;
    m_cumulativeEnergySquared /= m_totalSamplesTaken;
    m_energy    = m_cumulativeEnergy;
    m_variance  = m_cumulativeEnergySquared - m_cumulativeEnergy * m_cumulativeEnergy;
}

double Sampler::getEnergyAverage(int iteration) {
    return m_cumulativeEnergy / m_totalSamplesTaken;
}

double Sampler::getVariance(int iteration) {
    const double E  = m_cumulativeEnergy        / m_totalSamplesTaken;
    const double E2 = m_cumulativeEnergySquared / m_totalSamplesTaken;
    return E2 - E*E;
}

double Sampler::getEnergyBlockAverage(int iteration) {
    const double blockSize  = iteration - m_currentBlockStart;
    const double Eb         = m_cumulativeBlockEnergy        / (blockSize);
    const double Eb2        = m_cumulativeBlockEnergySquared / (blockSize);

    m_blockVariance                 = Eb2 - Eb*Eb;
    m_cumulativeBlockEnergy         = 0;
    m_cumulativeBlockEnergySquared  = 0;
    m_currentBlockStart             = iteration;
    return Eb;
}

double Sampler::getVarianceBlock(int) {
    const double Vb = m_blockVariance;
    m_blockVariance = 0;
    return Vb;
}

Sampler::Sampler(System* system) {
    m_system = system;
}

void Sampler::setup() {
    m_currentBlockStart = -1;
    m_hamiltonian = m_system->getHamiltonian();
}
