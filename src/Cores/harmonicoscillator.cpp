#include "harmonicoscillator.h"
#include "system.h"
#include "metropolis.h"
#include "electron.h"
#include <random>

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
    std::normal_distribution<double> normalDistributionX{m_position(0), m_size};
    std::normal_distribution<double> normalDistributionY{m_position(1), m_size};
    std::normal_distribution<double> normalDistributionZ{m_position(2), m_size};
    std::mt19937& randomGenerator = m_system->getMetropolis()->m_randomGenerator;
    for (int i = 0; i < m_numberOfElectrons; i++) {
        m_electrons.push_back(new Electron(vec {normalDistributionX(randomGenerator),
                                                normalDistributionY(randomGenerator),
                                                normalDistributionZ(randomGenerator)},
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


