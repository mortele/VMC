#pragma once

class Metropolis {
    friend class System;

private:
    int                 m_numberOfElectrons;
    int                 m_numberOfDimensions;
    double              m_stepLength = 1e-3;
    double              m_stepLengthHalf;
    class System*       m_system;
    class Sampler*      m_sampler;
    class WaveFunction* m_waveFunction;

    bool step();
    void setup();

public:
    Metropolis(class System* system);
    void runSteps(int steps);
};

