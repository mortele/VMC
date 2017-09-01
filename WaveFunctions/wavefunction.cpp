#include "wavefunction.h"
#include "system.h"

using arma::vec;
using arma::mat;
using arma::zeros;


WaveFunction::WaveFunction(System* system) {
    m_system = system;
    m_2stepLengthInverse = 1.0 / (2.0 * m_stepLength);
    m_stepLengthSquared = m_stepLength * m_stepLength;
    m_stepLengthSquaredInverse = 1.0 / m_stepLengthSquared;
    m_useNumericalDerivatives = false;
}

void WaveFunction::setup() {
    m_numberOfDimensions        = m_system->getNumberOfDimensions();
    m_numberOfElectrons         = m_system->getNumberOfElectrons();
    m_numberOfSpinUpElectrons   = m_system->getNumberOfSpinUpElectrons();
    m_numberOfSpinDownElectrons = m_system->getNumberOfSpinDownElectrons();

}

void WaveFunction::updateOldWaveFunctionValue() {
    m_oldValueSquared = m_currentValueSquared;
}

void WaveFunction::evaluateWaveFunctionInitial() {
    m_currentValueSquared = evaluateWaveFunctionSquared();
    m_oldValueSquared     = m_currentValueSquared;
}

double WaveFunction::computeWaveFunctionRatio(int) {
    m_currentValueSquared = evaluateWaveFunctionSquared();
    return m_currentValueSquared / m_oldValueSquared;
}

void WaveFunction::setStepLength(double stepLength) {
    m_stepLength                = stepLength;
    m_2stepLengthInverse        = 1.0 / (2.0 * m_stepLength);
    m_stepLengthSquared         = m_stepLength * m_stepLength;
    m_stepLengthSquaredInverse  = 1.0 / m_stepLengthSquared;
}

double WaveFunction::evaluateWaveFunctionSquared() {
    const double current = evaluateWaveFunction();
    m_currentValueSquared = current * current;
    return m_currentValueSquared;
}

double WaveFunction::evaluateLaplacian() {
    double laplacian = 0;
    double waveFunction = evaluateWaveFunction();

    for (int electron = 0; electron < m_numberOfElectrons; electron++) {
        for (int dimension = 0; dimension < m_numberOfDimensions; dimension++) {

            // Psi(x+dx)
            m_system->adjustPositionOfElectron(electron, dimension, m_stepLength);
            double waveFunctionPlus = evaluateWaveFunction();
            // Psi(x-dx)
            m_system->adjustPositionOfElectron(electron, dimension, -2*m_stepLength);
            double waveFunctionMinus = evaluateWaveFunction();
            // Adjust back to original position.
            m_system->adjustPositionOfElectron(electron, dimension, m_stepLength);

            // d^2/dx^2 Psi(x) = [ Psi(x+dx) - 2 Psi(x) + Psi(x-dx) ] / dx^2
            laplacian += waveFunctionPlus - 2* waveFunction + waveFunctionMinus;
        }
    }
    return laplacian * m_stepLengthSquaredInverse / waveFunction;
}

mat WaveFunction::evaluateGradient() {
    mat gradient = zeros<mat>(m_numberOfElectrons, m_numberOfDimensions);
    double waveFunction = evaluateWaveFunction();

    for (int electron = 0; electron < m_numberOfElectrons; electron++) {
        for (int dimension = 0; dimension < m_numberOfDimensions; dimension++) {

            // Psi(x+dx)
            m_system->adjustPositionOfElectron(electron, dimension, m_stepLength);
            double waveFunctionPlus = evaluateWaveFunction();

            // Psi(x-dx)
            m_system->adjustPositionOfElectron(electron, dimension, -2*m_stepLength);
            double waveFunctionMinus = evaluateWaveFunction();

            // d/dx Psi(x) = [ Psi(x+dx) - Psi(x-dx) ] / 2 dx
            m_system->adjustPositionOfElectron(electron, dimension, m_stepLength);
            gradient(electron, dimension) = waveFunctionPlus - waveFunctionMinus;
        }
    }
    gradient *= m_2stepLengthInverse / waveFunction;
    return gradient;
}
