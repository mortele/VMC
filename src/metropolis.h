#pragma once
#include <random>

class Metropolis {
    friend class System;
    friend class Atom;
    friend class HarmonicOscillator;

private:
    int                 m_i;
    int                 m_numberOfElectrons;
    int                 m_numberOfDimensions;
    int                 m_numberOfMetropolisSteps;
    bool                m_silent = false;
    bool                m_stepLengthSetManually = false;
    bool                m_importanceSampling = false;
    double              m_stepLength = 1.0;
    double              m_stepLengthHalf;
    double              m_dt;
    double              m_dtSqrt;
    class System*       m_system;
    class Sampler*      m_sampler;

protected:
    class WaveFunction* m_waveFunction;
    std::random_device  m_randomDevice;
    std::mt19937        m_randomGenerator{m_randomDevice()};

    bool step();
    void setup();
    void setStepLength(double stepLength);
    void printInitialInfo();
    void printIterationInfo(int iteration);
    void printFinalInfo();
    double computeGreensFunction();

public:
    Metropolis(class System* system);
    void setImportanceSampling(bool importanceSampling);
    double runSteps(int steps);
    double runStepsSilent(int steps);
    int  getStep() { return m_i; }
};

