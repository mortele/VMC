#pragma once
#include "sampler.h"

class Hamiltonian {
    friend class Sampler;
    friend class System;

protected:
    double              m_localEnergy;
    class System*       m_system;
    class WaveFunction* m_waveFunction;

    void   setup();
    double computeKineticEnergy();
    double computeElectronCorePotentialEnergy();
    double computeCoreCorePotentialEnergy();
    double computeElectronElectronPotentialEnergy();

public:
    Hamiltonian(class System* system);
    double computeLocalEnergy();
    double getLocalEnergy() { return m_localEnergy; }
};

