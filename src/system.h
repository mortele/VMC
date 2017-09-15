#pragma once
#include <vector>

class System {
private:
    int                             m_numberOfDimensions = 3;
    int                             m_numberOfElectrons;
    int                             m_numberOfSpinUpElectrons;
    int                             m_numberOfSpinDownElectrons;
    bool                            m_interactingElectrons = true;
    class Sampler*                  m_sampler;
    class Metropolis*               m_metropolis;
    class Hamiltonian*              m_hamiltonian;
    class WaveFunction*             m_waveFunction;
    class Orbital*                  m_orbital;
    std::vector<class Core*>        m_cores;
    std::vector<class Electron*>    m_electrons;
    std::vector<class Electron*>    m_spinUpElectrons;
    std::vector<class Electron*>    m_spinDownElectrons;

    void setup();
    void setupSpinElectronArrays();

public:
    System();
    double runMetropolis                (int steps);
    double runMetropolisSilent          (int steps);
    void adjustPositionOfElectron     (int electron, int dimension, double adjustment);
    void addCore                      (class Core* core);
    void setWaveFunction              (class WaveFunction* waveFunction);
    void setElectronInteraction       (bool interacting);
    void setStepLength                (double stepLength);
    void setImportanceSampling        (bool importanceSampling);
    void setOrbital                   (Orbital* orbital);
    int  getNumberOfDimensions        () { return m_numberOfDimensions;        }
    int  getNumberOfElectrons         () { return m_numberOfElectrons;         }
    int  getNumberOfSpinUpElectrons   () { return m_numberOfSpinUpElectrons;   }
    int  getNumberOfSpinDownElectrons () { return m_numberOfSpinDownElectrons; }
    class Sampler*               getSampler()           { return m_sampler;            }
    class Hamiltonian*           getHamiltonian()       { return m_hamiltonian;        }
    class WaveFunction*          getWaveFunction()      { return m_waveFunction;       }
    class Metropolis*            getMetropolis()        { return m_metropolis;         }
    std::vector<class Core*>     getCores()             { return m_cores;              }
    std::vector<class Electron*> getElectrons()         { return m_electrons;          }
    std::vector<class Electron*> getSpinUpElectrons()   { return m_spinUpElectrons;    }
    std::vector<class Electron*> getSpinDownElectrons() { return m_spinDownElectrons;  }
};

