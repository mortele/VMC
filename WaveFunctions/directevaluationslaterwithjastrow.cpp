#include "directevaluationslaterwithjastrow.h"
#include "system.h"
#include <cmath>
#include <armadillo>

using std::exp;
using arma::vec;
using arma::norm;


DirectEvaluationSlaterWithJastrow::DirectEvaluationSlaterWithJastrow(System* system,
                                                                     double alpha,
                                                                     double beta,
                                                                     int numberOfSpinUpElectrons,
                                                                     int numberOfSpinDownElectrons) :
        DirectEvaluationSlater(system,
                               alpha,
                               numberOfSpinUpElectrons,
                               numberOfSpinDownElectrons) {
    m_beta = beta;
    m_useNumericalDerivatives = true;

}

void DirectEvaluationSlaterWithJastrow::evaluateWaveFunctionInitial() {
    DirectEvaluationSlater::evaluateWaveFunctionInitial();
    const double jastrowFactor = computeJastrowFactor();
    m_currentValueSquared = m_currentValueSquared * jastrowFactor * jastrowFactor;
}

DirectEvaluationSlaterWithJastrow::evaluateWaveFunction() {
    const double currentWaveFunction = DirectEvaluationSlater::evaluateWaveFunction() *
                                       computeJastrowFactor();
    m_currentValueSquared = currentWaveFunction * currentWaveFunction;
    return currentWaveFunction;
}

double DirectEvaluationSlaterWithJastrow::computeJastrowFactor() {
    double jastrowLogarithm = 0;
    for (unsigned int electron1 = 0; electron1 < m_system->getElectrons().size(); electron1++) {
        for (unsigned int electron2 = electron1+1; electron2 < m_system->getElectrons().size(); electron2++) {
            const double r12 = norm(m_system->getElectrons().at(electron1)->getPosition() -
                                  m_system->getElectrons().at(electron2)->getPosition());
            const int spin1 = m_system->getElectrons().at(electron1)->getSpin();
            const int spin2 = m_system->getElectrons().at(electron2)->getSpin();
            const double spinFactor = (spin1 == spin2) ? 0.25 : 0.5;

            jastrowLogarithm += spinFactor * r12 / (1 + r12 * m_beta);
        }
    }
    return exp(jastrowLogarithm);
}
