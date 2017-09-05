#pragma once
#include "WaveFunctions/Orbitals/orbital.h"
#include <string>
#include <vector>
#include <armadillo>


class GaussianOrbital : public Orbital {
private:
    int         m_basisSize;
    arma::mat   m_spinUpCoefficients;
    arma::mat   m_spinDownCoefficients;
    class HartreeFockBasisParser* m_parser;
    std::vector<class ContractedGaussian*> m_basis;


public:
    GaussianOrbital(std::string basisFileName);
    double evaluate(double x, double y, double z, int index);
    double computeDerivative(double x, double y, double z, int index, int dimension);
    double computeLaplacian(double x, double y, double z, int index);
};
