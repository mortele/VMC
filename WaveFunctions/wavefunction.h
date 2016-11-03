#pragma once
#include <armadillo>

class WaveFunction {
protected:
    int             m_numberOfElectrons;
    int             m_numberOfDimensions;
    double          m_stepLength = 1e-4;
    double          m_stepLengthSquared;
    double          m_2stepLengthInverse;
    double          m_stepLengthSquaredInverse;
    class System*   m_system;

public:
    WaveFunction(class System* system);

    void setup();

    void setStepLength(double stepLength);

    virtual void        updateWaveFunction   () {}
    virtual double      evaluateWaveFunction () = 0;
    virtual double      evaluateLaplacian    ();
    virtual arma::mat   evaluateGradient     ();
};

