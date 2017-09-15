#pragma once
#include "WaveFunctions/heliumwavefunction.h"
#include <armadillo>

class HeliumWithJastrow : public HeliumWaveFunction {
private:
    double m_beta = 1.0;

public:
    HeliumWithJastrow(class System* system,
                      double        alpha,
                      double        beta,
                      bool          useNumericalDerivatives);
    double evaluateWaveFunction();
    void evaluateWaveFunctionInitial();
    double evaluateLaplacian();
};
