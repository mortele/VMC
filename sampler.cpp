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

    /*
     * http://journals.aps.org/prb/pdf/10.1103/PhysRevB.42.3503
     *
     * Should think about using :
     *
     *        1   --
     * <E> = ---  >  [ 1 - P(R -> R') ] E (R) + P(R -> R') E (R')
     *        N   --                     L                  L
     *
     * Same average as ordinary estimator, but smaller variance usually.
     */

    if (acceptedStep) {
        m_currentEnergy             = m_hamiltonian->computeLocalEnergy();
        m_currentKineticEnergy      = m_hamiltonian->getKineticEnergy();
        m_currentPotentialEnergy    = m_hamiltonian->getPotentialEnergy();
    } else {
        m_currentEnergy             = m_hamiltonian->getLocalEnergy();
        m_currentKineticEnergy      = m_hamiltonian->getKineticEnergy();
        m_currentPotentialEnergy    = m_hamiltonian->getPotentialEnergy();
    }

    m_currentEnergySquared            = m_currentEnergy * m_currentEnergy;
    m_cumulativeEnergy               += m_currentEnergy;
    m_cumulativeEnergySquared        += m_currentEnergySquared;
    m_blockCumulativeEnergy          += m_currentEnergy;
    m_blockCumulativeEnergySquared   += m_currentEnergySquared;
    m_blockCumulativeKineticEnergy   += m_currentKineticEnergy;
    m_blockCumulativePotentialEnergy += m_currentPotentialEnergy;
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
    m_variance                 = (E2-E*E) / m_totalSamplesTaken;
    m_numberOfMetropolisSteps  = m_totalSamplesTaken;
}

void Sampler::computeBlockAverages() {
    if (m_firstBlock) {
        m_firstBlock = false;
        m_currentEnergySquared     = 0;
        m_cumulativeEnergy         = 0;
        m_cumulativeEnergySquared  = 0;
        m_totalAccepted            = 0;
        m_totalSamplesTaken        = 0;
    }
    const double N        = m_blockSamplesTaken;
    const double E        = m_blockCumulativeEnergy             / N;
    const double E2       = m_blockCumulativeEnergySquared      / N;
    const double T        = m_blockCumulativeKineticEnergy      / N;
    const double V        = m_blockCumulativePotentialEnergy    / N;
    m_blockEnergy         = E;
    m_blockVariance       = (E2 - E*E) / m_blockSamplesTaken;

    m_blockAcceptanceRate = m_blockAccepted / N;
    m_blockVirialRatio    = T / V;

    m_blockSamplesTaken                 = 0;
    m_blockAccepted                     = 0;
    m_blockCumulativeEnergy             = 0;
    m_blockCumulativeEnergySquared      = 0;
    m_blockCumulativeKineticEnergy      = 0;
    m_blockCumulativePotentialEnergy    = 0;
}

Sampler::Sampler(System* system) {
    m_system = system;
}

void Sampler::setup() {
    m_hamiltonian = m_system->getHamiltonian();
}
