#pragma once
#include "WaveFunctions/directevaluationslater.h"

class DirectEvaluationSlaterWithJastrow : public DirectEvaluationSlater {
private:
    double m_beta;

    double computeJastrowFactor();

public:
    DirectEvaluationSlaterWithJastrow(class System* system,
                                      double        alpha,
                                      double        beta,
                                      int           numberOfSpinUpElectrons,
                                      int           numberOfSpinDownElectrons);
    double evaluateWaveFunction();
    void evaluateWaveFunctionInitial();
};
