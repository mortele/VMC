#pragma once
#include <armadillo>

class WaveFunction {
    friend class System;
    friend class Metropolis;

protected:
    int             m_numberOfElectrons;
    int             m_numberOfSpinUpElectrons;
    int             m_numberOfSpinDownElectrons;
    int             m_numberOfDimensions;
    int             m_changedElectron;
    int             m_changedDimension;
    bool            m_orbitalSet = false;
    bool            m_useNumericalDerivatives;
    bool            m_containsJastrow = false;
    double          m_currentValueSquared;
    double          m_oldValueSquared;
    double          m_stepLength = 1e-4;
    double          m_stepLengthSquared;
    double          m_2stepLengthInverse;
    double          m_stepLengthSquaredInverse;
    double          m_betaDerivative;
    arma::mat       m_quantumForce;
    arma::mat       m_quantumForceOld;
    arma::mat       m_electronPositions;
    arma::mat       m_electronPositionsOld;

    class System*   m_system;
    class Orbital*  m_orbital;

    virtual void    setup();
    void            updateOldWaveFunctionValue();
    virtual void    evaluateWaveFunctionInitial();
    virtual double  getQuantumForce   (int,int) { return nan(""); }
    virtual double  getQuantumForceOld(int,int) { return nan(""); }
    virtual double  computeWaveFunctionRatio(int changedElectronIndex);
    double          getPosition(int,int);
    double          getPositionOld(int,int);

public:
    WaveFunction(class System* system);
    void   setStepLength(double stepLength);
    void   setOrbital(class Orbital* orbital);
    double evaluateWaveFunctionSquared();
    virtual void setBeta(double) {}
    virtual double  computeBetaDerivative();

    bool   containsJastrow()   { return m_containsJastrow; }
    double getBetaDerivative() { return m_betaDerivative; }

    virtual double      evaluateWaveFunction () { return nan(""); }
    virtual void        passProposedChangeToWaveFunction(int , int ) {}
    virtual void        updateWaveFunctionAfterAcceptedStep() {}
    virtual void        updateWaveFunctionAfterRejectedStep() {}
    virtual double      evaluateLaplacian    ();
    virtual arma::mat   evaluateGradient     ();
};


inline double WaveFunction::getPosition(int i, int j) {
    return m_electronPositions(i,j);
}
inline double WaveFunction::getPositionOld(int i, int j) {
    return m_electronPositionsOld(i,j);
}
