#pragma once
#include "WaveFunctions/wavefunction.h"
#include <armadillo>

class HeliumWaveFunction : public WaveFunction {
private:
    bool    m_useNumericalDerivatives = false;
    double  m_alpha;

public:
    HeliumWaveFunction(class System* system,
                       double        alpha = 1.0,
                       bool          useNumericalDerivatives = false);
    void setup();
    double evaluateWaveFunction();
    double evaluateLaplacian();
    arma::mat evaluateGradient();

    void setAlpha(double alpha) { m_alpha = alpha; }
};

