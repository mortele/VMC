#include "gaussianslater.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/Orbitals/contractedgaussian.h"
#include "WaveFunctions/Orbitals/primitivegaussian.h"
#include "Cores/atom.h"
#include "system.h"
#include "electron.h"
#include <cmath>

using std::cout;
using std::endl;
using std::exp;
using arma::zeros;
using arma::mat;
using arma::det;
using arma::norm;


GaussianSlater::GaussianSlater(System* system, HartreeFockBasisParser* parser) :
        WaveFunction(system) {
    m_beta                      = 0.4;
    m_basis                     = parser->getBasis();
    m_basisSize                 = m_basis.size();
    m_numberOfDimensions        = 3;
    m_numberOfElectrons         = parser->getNumberOfElectrons();
    m_numberOfSpinUpElectrons   = parser->getNumberOfSpinUpElectrons();
    m_numberOfSpinDownElectrons = parser->getNumberOfSpinDownElectrons();
    m_spinUpCoefficients        = parser->getSpinUpCoefficients();
    m_spinDownCoefficients      = parser->getSpinDownCoefficients();

    for (int atom = 0; atom < parser->getNumberOfAtoms(); atom++) {
        int       Z         = parser->getAtomCharges().at(atom);
        arma::vec position  = parser->getAtomPosition().at(atom);
        m_system->addCore(new Atom(system, position, Z));
    }
}

double GaussianSlater::evaluateWaveFunction() {
    for (int electron = 0; electron < m_numberOfSpinUpElectrons; electron++) {
        for (int orbital = 0; orbital < m_numberOfSpinUpElectrons; orbital++) {
            m_spinUpDeterminant(electron, orbital) = computeSpinUpOrbital(electron, orbital);
        }
    }
    for (int electron = 0; electron < m_numberOfSpinDownElectrons; electron++) {
        for (int orbital = 0; orbital < m_numberOfSpinDownElectrons; orbital++) {
            m_spinDownDeterminant(electron, orbital) = computeSpinDownOrbital(electron, orbital);
        }
    }
    double logJastrow = 0;
    for (int electron1 = 0; electron1 < m_numberOfElectrons; electron1++) {
        for (int electron2 = electron1 + 1; electron2 < m_numberOfElectrons; electron2++) {
            const int       spin1   = m_system->getElectrons().at(electron1)->getSpin();
            const int       spin2   = m_system->getElectrons().at(electron2)->getSpin();
            const double    a       = 0.5 - 0.25 * (spin1==spin2);
            const double    r12     = norm(m_system->getElectrons().at(electron1)->getPosition() -
                                           m_system->getElectrons().at(electron2)->getPosition());
            logJastrow += a * r12 / (1 + m_beta * r12);
        }
    }
    //cout << logJastrow << " " << exp(logJastrow) << endl;
    return det(m_spinUpDeterminant) * det(m_spinDownDeterminant) * exp(logJastrow);
}

void GaussianSlater::evaluateWaveFunctionInitial() {
    m_spinUpDeterminant   = zeros<mat>(m_numberOfSpinUpElectrons, m_numberOfSpinUpElectrons);
    m_spinDownDeterminant = zeros<mat>(m_numberOfSpinDownElectrons, m_numberOfSpinDownElectrons);
    m_currentValueSquared = evaluateWaveFunction();
    m_currentValueSquared = m_currentValueSquared * m_currentValueSquared;
}

double GaussianSlater::evaluateLaplacian() {

}

double GaussianSlater::computeSpinUpOrbital(int electron, int orbital) {
    double result = 0;
    const double x = m_system->getElectrons().at(electron)->getPosition()(0);
    const double y = m_system->getElectrons().at(electron)->getPosition()(1);
    const double z = m_system->getElectrons().at(electron)->getPosition()(2);

    for (int basisFunction = 0; basisFunction < m_basisSize; basisFunction++) {
        result += m_spinUpCoefficients(basisFunction, orbital) *
                  (*m_basis.at(basisFunction))(x,y,z);
    }
    return result;
}

double GaussianSlater::computeSpinDownOrbital(int electron, int orbital) {
    double result = 0;
    const double x = m_system->getElectrons().at(electron)->getPosition()(0);
    const double y = m_system->getElectrons().at(electron)->getPosition()(1);
    const double z = m_system->getElectrons().at(electron)->getPosition()(2);

    for (int basisFunction = 0; basisFunction < m_basisSize; basisFunction++) {
        result += m_spinDownCoefficients(basisFunction, orbital) *
                  (*m_basis.at(basisFunction))(x,y,z);
    }
    return result;
}
