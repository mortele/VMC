#pragma once

class Sampler {
    friend class Metropolis;
    friend class System;

private:
    bool     m_firstBlock                        = true;
    int      m_numberOfMetropolisSteps           = 0;
    double   m_currentEnergy                     = 0;
    double   m_currentEnergySquared              = 0;
    double   m_currentKineticEnergy              = 0;
    double   m_currentPotentialEnergy            = 0;
    double   m_energy                            = 0;
    double   m_variance                          = 0;

    double   m_totalAccepted                     = 0;
    double   m_totalSamplesTaken                 = 0;
    double   m_cumulativeEnergy                  = 0;
    double   m_cumulativeEnergySquared           = 0;
    double   m_acceptanceRate                    = 0;

    double   m_blockEnergy                       = 0;
    double   m_blockVariance                     = 0;
    double   m_blockAccepted                     = 0;
    double   m_blockSamplesTaken                 = 0;
    double   m_blockCumulativeEnergy             = 0;
    double   m_blockCumulativeEnergySquared      = 0;
    double   m_blockAcceptanceRate               = 0;
    double   m_blockCumulativeKineticEnergy      = 0;
    double   m_blockCumulativePotentialEnergy    = 0;
    double   m_blockVirialRatio                  = 0;

    double   m_blockingEnergy                    = 0;
    double   m_blockingEnergy2                   = 0;
    double   m_blockingVariance                  = 0;
    double   m_blockingStandardDeviation         = 0;
    int      m_numberOfBlocks                    = 0;

    class System*       m_system      = nullptr;
    class Hamiltonian*  m_hamiltonian = nullptr;

    void setup();
    void sample(bool acceptedStep);
    void computeAverages();
    void computeBlockAverages();

    double getBlockEnergy()         { return m_blockEnergy;         }
    double getBlockVariance()       { return m_blockVariance;       }
    double getBlockAcceptanceRate() { return m_blockAcceptanceRate; }
    double getBlockVirialRatio()    { return m_blockVirialRatio;    }

public:
    Sampler(class System* system);

    double getEnergy()              { return m_energy;                      }
    double getVarianceBlocking()    { return m_blockingVariance;            }
    double getVariance()            { return m_variance;                    }
    double getStandardDeviation()   { return m_blockingStandardDeviation;   }
    double getAcceptanceRate()      { return m_acceptanceRate;              }
};
