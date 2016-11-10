#include "directevaluationslater.h"
#include "system.h"
#include "electron.h"
#include <iostream>
#include <cmath>
#include <cassert>

using std::cout;
using std::endl;
using std::exp;
using arma::norm;
using arma::zeros;
using arma::mat;
using arma::det;

DirectEvaluationSlater::DirectEvaluationSlater(System*  system,
                                               double   alpha,
                                               int      numberOfSpinUpElectrons,
                                               int      numberOfSpinDownElectrons,
                                               bool     useNumericalDerivatives) :
        WaveFunction(system) {
    m_numberOfSpinUpElectrons   = numberOfSpinUpElectrons;
    m_numberOfSpinDownElectrons = numberOfSpinDownElectrons;
    m_numberOfSpatialOrbitals   = numberOfSpinUpElectrons;
    m_useNumericalDerivatives   = useNumericalDerivatives;
    m_alpha                     = alpha;

    assert(numberOfSpinUpElectrons == numberOfSpinDownElectrons);
}

double DirectEvaluationSlater::evaluateWaveFunction() {
    const double determinantProduct = m_spinUpDeterminant * m_spinDownDeterminant;
}

void DirectEvaluationSlater::evaluateWaveFunctionInitial() {
    m_slaterSpinUp   = zeros<mat>(m_numberOfSpinUpElectrons,   m_numberOfSpatialOrbitals);
    m_slaterSpinDown = zeros<mat>(m_numberOfSpinDownElectrons, m_numberOfSpatialOrbitals);

    for (int electron = 0; electron < m_numberOfSpinUpElectrons; electron++) {
        for (int orbital = 0; orbital  < m_numberOfSpatialOrbitals; orbital ++) {
            m_slaterSpinUp  (electron, orbital) = evaluateOrbital(orbital, electron, 1);
            m_slaterSpinDown(electron, orbital) = evaluateOrbital(orbital, electron, 0);
        }
    }
    m_spinUpDeterminant             = det(m_slaterSpinUp);
    m_spinDownDeterminant           = det(m_slaterSpinDown);
    const double determinantProduct = m_spinUpDeterminant * m_spinDownDeterminant;
    m_currentValueSquared           = determinantProduct * determinantProduct;
}

double DirectEvaluationSlater::computeWaveFunctionRatio(int changedElectronIndex) {
    int spin = m_system->getElectrons().at(changedElectronIndex)->getSpin();
    if (spin == 1) {
        updateSpinUpSlater();
    } else if (spin == 0) {
        updateSpinDownSlater();
    } else {
        std::cout << "Something is horribly wrong in DirectEvaluationSlater::computeWaveFunctionRatio(int)." << std::endl;
    }
    const double waveFunction = evaluateWaveFunction();
    m_currentValueSquared = waveFunction * waveFunction;
    return m_currentValueSquared / m_oldValueSquared;
}

double DirectEvaluationSlater::evaluateLaplacian() {

}

double DirectEvaluationSlater::evaluateOrbital(int orbital, int electron, bool up) {
    if (orbital == 0) {
        return s1(electron, up);
    } else if (orbital == 1) {
        return s2(electron, up);
    } else if (orbital == 2) {
        return p2x(electron, up);
    } else if (orbital == 3) {
        return p2y(electron, up);
    } else if (orbital == 4) {
        return p2z(electron, up);
    } else {
        cout << "Unknown orbital #" << orbital << endl;
        return NAN;
    }
}

void DirectEvaluationSlater::updateSpinUpSlater() {
    for (int electron = 0; electron < m_numberOfSpinUpElectrons; electron++) {
        for (int orbital = 0; orbital  < m_numberOfSpatialOrbitals; orbital ++) {
            m_slaterSpinUp  (electron, orbital) = evaluateOrbital(orbital, electron, 1);
        }
    }
    m_spinUpDeterminant = det(m_slaterSpinUp);
}

void DirectEvaluationSlater::updateSpinDownSlater() {
    for (int electron = 0; electron < m_numberOfSpinDownElectrons; electron++) {
        for (int orbital = 0; orbital  < m_numberOfSpatialOrbitals; orbital ++) {
            m_slaterSpinDown(electron, orbital) = evaluateOrbital(orbital, electron, 0);
        }
    }
    m_spinDownDeterminant = det(m_slaterSpinDown);
}

double DirectEvaluationSlater::s1(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    return exp(-m_alpha * r);
}

double DirectEvaluationSlater::s2(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    return (1 - 0.5 * m_alpha * r) * exp(-0.5 * m_alpha * r);
}

double DirectEvaluationSlater::p2x(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    const double x = electrons.at(electron)->getPosition()(0);
    return x * exp(-m_alpha * r);
}

double DirectEvaluationSlater::p2y(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    const double y = electrons.at(electron)->getPosition()(1);
    return y * exp(-m_alpha * r);
}

double DirectEvaluationSlater::p2z(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    const double z = electrons.at(electron)->getPosition()(2);
    return z * exp(-m_alpha * r);

}

double DirectEvaluationSlater::s1Laplacian(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    return (m_alpha * m_alpha - 2 * m_alpha / r) * exp(-m_alpha * r);
}

double DirectEvaluationSlater::s2Laplacian(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    return (  5./4 * m_alpha * m_alpha
            - 2 * m_alpha / r
            - m_alpha * m_alpha * m_alpha / 8.) * exp(-0.5 * m_alpha * r);
}

double DirectEvaluationSlater::p2xLaplacian(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    const double x = electrons.at(electron)->getPosition()(0);
    return (m_alpha * x * (m_alpha * r - 8.)) / (4. * r) * exp(-0.5 * m_alpha * r);
}

double DirectEvaluationSlater::p2yLaplacian(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    const double y = electrons.at(electron)->getPosition()(1);
    return (m_alpha * y * (m_alpha * r - 8.)) / (4. * r) * exp(-0.5 * m_alpha * r);
}

double DirectEvaluationSlater::p2zLaplacian(int electron, bool up) {
    std::vector<Electron*> electrons = (up == 1) ? m_system->getSpinUpElectrons() :
                                                   m_system->getSpinDownElectrons();
    const double r = norm(electrons.at(electron)->getPosition());
    const double z = electrons.at(electron)->getPosition()(2);
    return (m_alpha * z * (m_alpha * r - 8.)) / (4. * r) * exp(-0.5 * m_alpha * r);
}





