#include "gaussianslater.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/Orbitals/contractedgaussian.h"
#include "WaveFunctions/Orbitals/primitivegaussian.h"
#include "Cores/atom.h"
#include "system.h"
#include "electron.h"

using arma::zeros;
using arma::mat;
using arma::det;


GaussianSlater::GaussianSlater(System* system, HartreeFockBasisParser* parser) :
        WaveFunction(system) {
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
    return det(m_spinUpDeterminant) * det(m_spinDownDeterminant);
}

void GaussianSlater::evaluateWaveFunctionInitial() {
    m_spinUpDeterminant   = zeros<mat>(m_numberOfSpinUpElectrons, m_numberOfSpinUpElectrons);
    m_spinDownDeterminant = zeros<mat>(m_numberOfSpinDownElectrons, m_numberOfSpinDownElectrons);
    m_oldValueSquared = evaluateWaveFunction();
    m_oldValueSquared = m_oldValueSquared * m_oldValueSquared;
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
