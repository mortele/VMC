#pragma once

class Metropolis {
private:
    class System*       m_system;
    class Sampler*      m_sampler;
    class WaveFunction* m_waveFunction;

    bool step();

public:
    Metropolis(class System* system);
    void setup();
    void runSteps(int steps);
};

