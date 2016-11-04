#include "sampler.h"
#include "system.h"
#include "hamiltonian.h"

void Sampler::sample(bool acceptedStep) {
    if (acceptedStep) {
        m_currentEnergy = m_hamiltonian->computeLocalEnergy();
    } else {
        m_currentEnergy = m_hamiltonian->getLocalEnergy();
    }
    m_currentEnergySquared     = m_currentEnergy * m_currentEnergy;
    m_cumulativeEnergy        += m_currentEnergy;
    m_cumulativeEnergySquared += m_currentEnergySquared;
}

Sampler::Sampler(System* system) {
    m_system = system;
}

void Sampler::setup() {
    m_hamiltonian = m_system->getHamiltonian();
}
