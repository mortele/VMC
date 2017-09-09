#pragma once

class Metropolis {
    friend class System;

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
    class WaveFunction* m_waveFunction;

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
    void runSteps(int steps);
    void runStepsSilent(int steps);
    int  getStep() { return m_i; }
};

