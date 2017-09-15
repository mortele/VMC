#pragma once
#include "sampler.h"

class Hamiltonian {
    friend class Sampler;
    friend class System;

protected:
    bool                m_interactingElectrons = true;
    double              m_localEnergy;
    double              m_kineticEnergy;
    double              m_potentialEnergy;
    class System*       m_system;
    class WaveFunction* m_waveFunction;

    void   setup();
    void   setElectronInteraction(bool interaction);
    double computeKineticEnergy();
    double computeElectronCorePotentialEnergy();
    double computeCoreCorePotentialEnergy();
    double computeElectronElectronPotentialEnergy();

public:
    Hamiltonian(class System* system);
    double computeLocalEnergy();
    double getLocalEnergy()     { return m_localEnergy; }
    double getKineticEnergy()   { return m_kineticEnergy; }
    double getPotentialEnergy() { return m_potentialEnergy; }
};

