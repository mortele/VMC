#pragma once
#include "wavefunction.h"
#include "WaveFunctions/Orbitals/hydrogenorbital.h"
#include "electron.h"
#include <armadillo>
#include <vector>


class System;


class SlaterWithJastrow : public WaveFunction {
private:
    bool        m_jastrow                   = true;
    int         m_spinChanged               = -1;
    double      m_alpha                     = 1;
    double      m_beta                      = 1;
    double      m_energyCrossTerm;
    double      m_jastrowLaplacian;
    double      m_slaterLaplacian;
    double      m_slaterLaplacianUp;
    double      m_slaterLaplacianDown;
    double      m_laplacian;
    double      m_Rsd;
    double      m_Rc;
    arma::mat   m_slaterUp;
    arma::mat   m_slaterDown;
    arma::mat   m_slaterGradientUp;
    arma::mat   m_slaterGradientDown;
    arma::mat   m_jastrowGradient;
    arma::mat   m_jastrowLaplacianTerms;
    arma::mat   m_positionsOld;
    arma::mat   m_correlationMatrix;
    arma::mat   m_correlationMatrixOld;

    HydrogenOrbital* m_orbital;

    void updateSlaterGradient(double Rsd, int electron);
    void updateJastrowGradient(int electron);
    void updateJastrowLaplacianTerms(int electron);
    void computeJastrowLaplacian();
    //void updateSlaterInverse();
    void computeSlaterRatio();
    void computeJastrowRatio();
    void computeSlaterLaplacian(int electron);
    void fillCorrelationMatrix();
    void updateCorrelationsMatrix();
    double computeInterEletronDistance(Electron*,Electron*);
    double computeJastrowFactor(Electron*,Electron*);
    double spinCoefficient(Electron*,Electron*);

protected:
    void computeQuantumForce();

public:
    SlaterWithJastrow(System* system, double alpha, double beta);
    void evaluateWaveFunctionInitial();
    void passProposedChangeToWaveFunction(int electronChanged, int dimensionChanged);
    void updateWaveFunctionAfterAcceptedStep();
    double computeWaveFunctionRatio(int changedElectronIndex);
    double evaluateLaplacian();
};


// Defined in header to make it possible for compiler optimization to inline it
// (hopefully... maybe...).
inline double SlaterWithJastrow::computeInterEletronDistance(Electron* a, Electron* b) {
    const double x = a->getPosition().at(0) - b->getPosition().at(0);
    const double y = a->getPosition().at(1) - b->getPosition().at(1);
    const double z = a->getPosition().at(2) - b->getPosition().at(2);
    return std::sqrt(x*x + y*y + z*z);
}
inline double SlaterWithJastrow::computeJastrowFactor(Electron* a,Electron* b) {
    const double aa     = spinCoefficient(a,b);
    const double rik    = computeInterEletronDistance(a,b);
    return aa * rik / (1 + m_beta * rik);
}
inline double SlaterWithJastrow::spinCoefficient(Electron* a, Electron* b) {
    return (a->getSpin() == b->getSpin() ? 0.25 : 0.5);
}









