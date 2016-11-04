#pragma once

class Sampler {
    friend class Metropolis;

private:
    double              m_currentEnergy;
    double              m_currentEnergySquared;
    double              m_cumulativeEnergy;
    double              m_cumulativeEnergySquared;
    class System*       m_system;
    class Hamiltonian*  m_hamiltonian;

    void sample(bool acceptedStep);

public:
    Sampler(class System* system);
    void setup();
};
