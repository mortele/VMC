#pragma once
#include "WaveFunctions/wavefunction.h"
#include <vector>
#include <armadillo>


class HydrogenWaveFunction : public WaveFunction {
private:
    double m_alpha;

public:
    HydrogenWaveFunction(class System* system,
                         double        alpha = 1.0,
                         bool          useNumericalDerivatives = false);
    void setup();
    double      evaluateWaveFunction();
    void evaluateWaveFunctionInitial();
    double      evaluateLaplacian();
    arma::mat   evaluateGradient();

    void setAlpha(double alpha) { m_alpha = alpha; }
};

