#include "harmonicoscillator.h"
#include "system.h"
#include "electron.h"
#include "RandomNumberGenerator/random.h"

using arma::vec;


HarmonicOscillator::HarmonicOscillator(System*      system,
                                       arma::vec    position,
                                       int          numberOfElectrons,
                                       double       omega) :
        Core(system, position) {
    m_generalizedCharge = omega;
    m_numberOfElectrons = numberOfElectrons;
    createElectrons();
}


void HarmonicOscillator::createElectrons() {
    for (int i = 0; i < m_numberOfElectrons; i++) {
        m_electrons.push_back(new Electron(vec {Random::nextGaussian(m_position(0), 1),
                                                Random::nextGaussian(m_position(1), 1),
                                                Random::nextGaussian(m_position(2), 1)},
                                           i % 2));
    }
}

double HarmonicOscillator::computeElectronCoreInteraction() {
    double interactionEnergy = 0;
    for (Electron* electron : m_system->getElectrons()) {
        arma::vec dr        = electron->getPosition() - m_position;
        double rSquared     = arma::dot(dr, dr);
        interactionEnergy  += 0.5 * m_generalizedCharge * rSquared;
    }
    return interactionEnergy;
}


