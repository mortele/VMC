#include "hydrogenwavefunction.h"
#include "system.h"
#include "electron.h"
#include <cmath>
#include <armadillo>
#include <cassert>

using arma::vec;
using arma::mat;
using arma::norm;
using arma::zeros;
using std::exp;

HydrogenWaveFunction::HydrogenWaveFunction(System* system,
                                           double  alpha,
                                           bool    useNumericalDerivatives) :
        WaveFunction(system) {
    m_alpha                     = alpha;
    m_useNumericalDerivatives   = useNumericalDerivatives;
}

void HydrogenWaveFunction::setup() {
    WaveFunction::setup();
    assert(m_numberOfElectrons == 1);
    assert(m_numberOfDimensions == 3);
}

double HydrogenWaveFunction::evaluateWaveFunction() {
    vec position = m_system->getElectrons().at(0)->getPosition();
    double r = norm(position);
    return exp(-m_alpha * r);
}

double HydrogenWaveFunction::evaluateLaplacian() {
    if (! m_useNumericalDerivatives) {
        vec     position = m_system->getElectrons().at(0)->getPosition();
        double  rInverse = 1.0 / norm(position);
        return m_alpha * m_alpha - (2.0 * m_alpha * rInverse);

    } else {
        return WaveFunction::evaluateLaplacian();
    }
}

arma::mat HydrogenWaveFunction::evaluateGradient() {
    if (! m_useNumericalDerivatives) {
        mat     gradient = zeros<mat>(m_numberOfElectrons, m_numberOfDimensions);
        vec     position = m_system->getElectrons().at(0)->getPosition();
        double  rInverse = 1.0 / norm(position);
        double  constant = - m_alpha * rInverse;
        for (int dimension = 0; dimension < m_numberOfDimensions; dimension++) {
            gradient(dimension) = constant * position(dimension);
        }
        return gradient;

    } else {
        return WaveFunction::evaluateGradient();
    }
}
