#pragma once
#include "sampler.h"

class Hamiltonian {
    friend class Sampler;

protected:
    double              m_localEnergy;
    class System*       m_system;
    class WaveFunction* m_waveFunction;

    double computeKineticEnergy();

public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy() = 0;
    double getLocalEnergy() { return m_localEnergy; }
};

