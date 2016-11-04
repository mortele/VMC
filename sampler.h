#pragma once

class Sampler {
    friend class Metropolis;
    friend class System;

private:
    int                 m_numberOfMetropolisSteps;
    int                 m_currentBlockStart;
    double              m_totalSamplesTaken;
    double              m_currentEnergy;
    double              m_currentEnergySquared;
    double              m_cumulativeEnergy;
    double              m_cumulativeEnergySquared;
    double              m_cumulativeBlockEnergy;
    double              m_cumulativeBlockEnergySquared;
    double              m_blockVariance;
    double              m_energy;
    double              m_variance;
    class System*       m_system;
    class Hamiltonian*  m_hamiltonian;

    void setup();
    void sample(bool acceptedStep);
    void computeAverages(int steps);
    double getEnergyAverage(int iteration);
    double getVariance(int iteration);
    double getEnergyBlockAverage(int iteration);
    double getVarianceBlock(int iteration);
    double getEnergy()   { return m_energy;   }
    double getVariance() { return m_variance; }

public:
    Sampler(class System* system);
};
