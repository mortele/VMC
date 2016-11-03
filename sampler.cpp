#include "sampler.h"
#include "system.h"

Sampler::Sampler(System* system) {
    m_system = system;
}

void Sampler::setup() {
    m_hamiltonian = m_system->getHamiltonian();
}
