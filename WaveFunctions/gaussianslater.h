#pragma once
#include "WaveFunctions/wavefunction.h"
#include <vector>
#include <armadillo>


class GaussianSlater : public WaveFunction {
private:
    int         m_basisSize;
    int         m_numberOfSpinUpElectrons;
    int         m_numberOfSpinDownElectrons;
    arma::mat   m_spinUpCoefficients;
    arma::mat   m_spinDownCoefficients;
    arma::mat   m_spinUpDeterminant;
    arma::mat   m_spinDownDeterminant;
    std::vector<class ContractedGaussian*> m_basis;

    double computeSpinUpOrbital  (int electron, int orbital);
    double computeSpinDownOrbital(int electron, int orbital);

public:
    GaussianSlater(class System* system, class HartreeFockBasisParser* parser);
    double evaluateWaveFunction();
    void evaluateWaveFunctionInitial();
};
