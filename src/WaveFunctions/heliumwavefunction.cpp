#include "heliumwavefunction.h"
#include "system.h"
#include "electron.h"
#include <cassert>

using arma::vec;
using arma::mat;
using arma::norm;
using arma::zeros;
using std::exp;

HeliumWaveFunction::HeliumWaveFunction(System* system,
                                       double  alpha,
                                       bool    useNumericalDerivatives) :
        WaveFunction(system) {
    m_useNumericalDerivatives = useNumericalDerivatives;
    m_alpha = alpha;
}

void HeliumWaveFunction::setup() {
    WaveFunction::setup();
    assert(m_numberOfElectrons  == 2);
    assert(m_numberOfDimensions == 3);
}

double HeliumWaveFunction::evaluateWaveFunction() {
    vec position1 = m_system->getElectrons().at(0)->getPosition();
    vec position2 = m_system->getElectrons().at(1)->getPosition();
    double r1 = norm(position1);
    double r2 = norm(position2);
    return exp(-m_alpha * (r1 + r2));
}

void HeliumWaveFunction::evaluateWaveFunctionInitial() {
    evaluateLaplacian();
}

double HeliumWaveFunction::evaluateLaplacian() {
    if (! m_useNumericalDerivatives) {
        vec position1 = m_system->getElectrons().at(0)->getPosition();
        vec position2 = m_system->getElectrons().at(1)->getPosition();
        const double r1Inverse = 1.0 / norm(position1);
        const double r2Inverse = 1.0 / norm(position2);
        return 2.0 * m_alpha * m_alpha - 2.0 * m_alpha * (r1Inverse + r2Inverse);
    } else {
        return WaveFunction::evaluateLaplacian();
    }
}

arma::mat HeliumWaveFunction::evaluateGradient() {
    if (! m_useNumericalDerivatives) {
        mat     gradient = zeros<mat>(m_numberOfElectrons, m_numberOfDimensions);
        for (int electron = 0; electron < m_numberOfElectrons; electron++) {
            vec     position = m_system->getElectrons().at(electron)->getPosition();
            double  rInverse = 1.0 / norm(position);
            double  constant = - m_alpha * rInverse;
            for (int dimension = 0; dimension < m_numberOfDimensions; dimension++) {
                gradient(electron, dimension) = constant * position(dimension);
            }
        }
        return gradient;
    } else {
        return WaveFunction::evaluateGradient();
    }
}






