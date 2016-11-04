#include "metropolis.h"
#include "system.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"

bool Metropolis::step() {
    return false;
}

Metropolis::Metropolis(System* system) {
    m_system = system;
}

void Metropolis::setup() {
    m_waveFunction  = m_system->getWaveFunction();
    m_sampler       = m_system->getSampler();
}

void Metropolis::runSteps(int steps) {
    for (int i = 0; i < steps; i++) {
        bool acceptedStep = step();
        m_sampler->sample(acceptedStep);
    }
}
