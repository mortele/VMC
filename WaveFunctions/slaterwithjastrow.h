#pragma once
#include "wavefunction.h"
#include "WaveFunctions/Orbitals/hydrogenorbital.h"
#include <armadillo>



class System;

class SlaterWithJastrow : public WaveFunction {
private:
    int         m_spinChanged = -1;
    arma::mat   m_slaterUp;
    arma::mat   m_slaterDown;
    arma::mat   m_slaterGradientUp;
    arma::mat   m_slaterGradientDown;
    arma::mat   m_jastrowGradient;

    HydrogenOrbital* m_orbital;

public:
    SlaterWithJastrow(System* system, double alpha, double beta);
    void evaluateWaveFunctionInitial();
    void updateWaveFunction(int electronChanged, int dimensionChanged);
};
