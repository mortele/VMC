#pragma once
#include <armadillo>

class WaveFunction {
    friend class System;
    friend class Metropolis;

protected:
    int             m_numberOfElectrons;
    int             m_numberOfDimensions;
    bool            m_useNumericalDerivatives;
    double          m_currentValueSquared;
    double          m_oldValueSquared;
    double          m_stepLength = 1e-4;
    double          m_stepLengthSquared;
    double          m_2stepLengthInverse;
    double          m_stepLengthSquaredInverse;
    class System*   m_system;

    void setup();
    void            updateOldWaveFunctionValue();
    virtual void    evaluateWaveFunctionInitial();
    virtual double  computeWaveFunctionRatio(int changedElectronIndex);

public:
    WaveFunction(class System* system);
    void   setStepLength(double stepLength);
    double evaluateWaveFunctionSquared ();

    virtual double      evaluateWaveFunction () = 0;
    virtual void        updateWaveFunction   () {}
    virtual double      evaluateLaplacian    ();
    virtual arma::mat   evaluateGradient     ();
};

