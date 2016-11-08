#pragma once
#include "WaveFunctions/wavefunction.h"
#include <armadillo>

class DirectEvaluationSlater : public WaveFunction {
private:
    int         m_numberOfSpinUpElectrons;
    int         m_numberOfSpinDownElectrons;
    double      m_alpha;
    arma::mat   m_slaterSpinUp;
    arma::mat   m_slaterSpinDown;

    double evaluateOrbital(int orbital, int electron);
    double s1(int electron);
    double s2(int electron);
    double p2x(int electron);
    double p2y(int electron);
    double p2z(int electron);

public:
    DirectEvaluationSlater(class System* system,
                           double        alpha,
                           int           numberOfSpinUpElectrons,
                           int           numberOfSpinDownElectrons,
                           bool          useNumericalDerivatives);
};

