#pragma once

class Optimizer {
    friend class System;

private:
    double              m_tollerance;
    double              m_cycles;
    double              m_maximumIterations;
    double              m_beta;
    bool                m_jastrow = false;

    class System*       m_system;
    class WaveFunction* m_waveFunction;
    class Metropolis*   m_metropolis;
    class Sampler*      m_sampler;

    void    printInitialInfo();
    void    printIterationInfo(int);
    void    setup();
    double  optimizeBeta(double beta, double tollerance, int maxIterations, int cycles);

public:
    Optimizer(class System* system);
};
