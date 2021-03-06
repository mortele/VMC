#pragma once
#include "WaveFunctions/wavefunction.h"
#include <armadillo>

class DirectEvaluationSlater : public WaveFunction {
protected:
    int         m_numberOfSpinUpElectrons;
    int         m_numberOfSpinDownElectrons;
    int         m_numberOfSpatialOrbitals;
    double      m_alpha;
    double      m_spinUpDeterminant;
    double      m_spinDownDeterminant;
    arma::mat   m_slaterSpinUp;
    arma::mat   m_slaterSpinDown;

    double evaluateOrbital(int orbital, int electron, bool up);
    void updateSpinUpSlater();
    void updateSpinDownSlater();

    double s1 (int electron, bool up);
    double s2 (int electron, bool up);
    double p2x(int electron, bool up);
    double p2y(int electron, bool up);
    double p2z(int electron, bool up);

    double s1Laplacian (int electron, bool up);
    double s2Laplacian (int electron, bool up);
    double p2xLaplacian(int electron, bool up);
    double p2yLaplacian(int electron, bool up);
    double p2zLaplacian(int electron, bool up);

    arma::mat computeSpinUpGradient();

public:
    DirectEvaluationSlater(class System* system,
                           double        alpha,
                           int           numberOfSpinUpElectrons,
                           int           numberOfSpinDownElectrons);
    double      evaluateWaveFunction();
    void        evaluateWaveFunctionInitial();
    double      computeWaveFunctionRatio(int changedElectronIndex);
    double      evaluateLaplacian();
};

