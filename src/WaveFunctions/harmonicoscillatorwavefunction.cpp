#include "harmonicoscillatorwavefunction.h"
#include "system.h"
#include "electron.h"
#include <armadillo>
#include <cmath>

HarmonicOscillatorWaveFunction::HarmonicOscillatorWaveFunction(System*  system,
                                                               double   alpha,
                                                               double   beta,
                                                               double   omega) :
        WaveFunction(system) {
    m_alpha = alpha;
    m_beta  = beta;
    m_omega = omega;
}

double HarmonicOscillatorWaveFunction::evaluateWaveFunction() {
    arma::vec position1 = m_system->getElectrons().at(0)->getPosition();
    arma::vec position2 = m_system->getElectrons().at(1)->getPosition();
    double r1Squared = arma::dot(position1, position1);
    double r2Squared = arma::dot(position2, position2);
    double r12       = arma::norm(position1 - position2);

    return std::exp(- 0.5 * m_omega * m_alpha * (r1Squared + r2Squared))
            * (m_beta == 0 ? 1.0 : std::exp(r12 / (2 + 2 * m_beta * r12)));
}

void HarmonicOscillatorWaveFunction::evaluateWaveFunctionInitial() {
    evaluateLaplacian();
}

double HarmonicOscillatorWaveFunction::evaluateLaplacian() {
    return WaveFunction::evaluateLaplacian();
}

