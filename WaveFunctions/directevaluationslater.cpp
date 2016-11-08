#include "directevaluationslater.h"
#include "system.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::exp;
using arma::norm;

DirectEvaluationSlater::DirectEvaluationSlater(System*  system,
                                               double   alpha,
                                               int      numberOfSpinUpElectrons,
                                               int      numberOfSpinDownElectrons,
                                               bool     useNumericalDerivatives) :
        WaveFunction(system) {
    m_numberOfSpinUpElectrons   = numberOfSpinUpElectrons;
    m_numberOfSpinDownElectrons = numberOfSpinDownElectrons;
    m_useNumericalDerivatives   = useNumericalDerivatives;
    m_alpha                     = alpha;

}
double DirectEvaluationSlater::evaluateWaveFunction() {

}

double DirectEvaluationSlater::evaluateOrbital(int orbital, int electron) {
    if (orbital == 0) {
        return s1(electron);
    } else if (orbital == 1) {
        return s2(electron);
    } else if (orbital == 2) {
        return p2x(electron);
    } else if (orbital == 3) {
        return p2y(electron);
    } else if (orbital == 4) {
        return p2z(electron);
    } else {
        cout << "Unknown orbital #" << orbital << endl;
        return nan;
    }
}

double DirectEvaluationSlater::s1(int electron) {
    const double r = norm(m_system->getElectrons().at(electron)->getPosition());
    return exp(-m_alpha * r);
}

double DirectEvaluationSlater::s2(int electron) {
    const double r = norm(m_system->getElectrons().at(electron)->getPosition());
    return (1 - 0.5 * m_alpha * r) * exp(-0.5 * m_alpha * r);
}

double DirectEvaluationSlater::p2x(int electron) {
    const double x = m_system->getElectrons().at(electron)->getPosition().at(0);
    const double r = norm(m_system->getElectrons().at(electron)->getPosition());
    return x * exp(-m_alpha * r);
}

double DirectEvaluationSlater::p2y(int electron) {
    const double y = m_system->getElectrons().at(electron)->getPosition().at(1);
    const double r = norm(m_system->getElectrons().at(electron)->getPosition());
    return y * exp(-m_alpha * r);
}

double DirectEvaluationSlater::p2z(int electron) {
    const double z = m_system->getElectrons().at(electron)->getPosition().at(2);
    const double r = norm(m_system->getElectrons().at(electron)->getPosition());
    return z * exp(-m_alpha * r);
}
