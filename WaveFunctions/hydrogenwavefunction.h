#pragma once
#include "WaveFunctions/wavefunction.h"
#include <vector>
#include <armadillo>


class HydrogenWaveFunction : public WaveFunction {
private:
    bool                            m_useNumericalDerivatives = true;
    double                          m_alpha;
    std::vector<class Electron*>    m_electrons;

public:
    HydrogenWaveFunction(class System* system, bool useNumericalDerivatives);
    void setup();
    double      evaluateWaveFunction();
    double      evaluateLaplacian();
    arma::mat   evaluateGradient();

    void setAlpha(double alpha) { m_alpha = alpha; }
};

