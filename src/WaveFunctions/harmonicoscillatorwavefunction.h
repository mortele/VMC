#pragma once
#include "WaveFunctions/wavefunction.h"


class HarmonicOscillatorWaveFunction : public WaveFunction {
private:
    double m_alpha;
    double m_beta;
    double m_omega;

public:
    HarmonicOscillatorWaveFunction(class System*    system,
                                   double           alpha,
                                   double           beta,
                                   double           omega);
    double evaluateWaveFunction();
    double evaluateLaplacian();
};
