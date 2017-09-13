#include "hamiltonian.h"
#include "WaveFunctions/hydrogenwavefunction.h"
#include "Cores/core.h"
#include "electron.h"
#include "system.h"
#include <armadillo>

using std::cout;
using std::endl;
using arma::vec;
using arma::norm;


void Hamiltonian::setup() {
    m_waveFunction = m_system->getWaveFunction();
}

void Hamiltonian::setElectronInteraction(bool interaction) {
    m_interactingElectrons = interaction;
}

double Hamiltonian::computeKineticEnergy() {
    return - 0.5 * m_waveFunction->evaluateLaplacian();
}

double Hamiltonian::computeElectronCorePotentialEnergy() {
    double electronCoreInteractionEnergy = 0;
    for (Core* core : m_system->getCores()) {
        electronCoreInteractionEnergy += core->computeElectronCoreInteraction();
    }
    return electronCoreInteractionEnergy;
}

double Hamiltonian::computeCoreCorePotentialEnergy() {
    double coreCoreInteractionEnergy = 0;
    for (Core* core : m_system->getCores()) {
        coreCoreInteractionEnergy += core->computeCoreCoreInteraction();
    }
    return 0.5 * coreCoreInteractionEnergy;
}

double Hamiltonian::computeElectronElectronPotentialEnergy() {
    if (m_interactingElectrons) {
        double electronElectronInteractionEnergy = 0;
        for (unsigned int i = 0; i < m_system->getElectrons().size(); i++) {
            for (unsigned int j = i+1; j < m_system->getElectrons().size(); j++) {
                //vec dr = m_system->getElectrons().at(i)->getPosition() -
                //         m_system->getElectrons().at(j)->getPosition();
                const double dx = m_system->getElectrons().at(i)->getPosition().at(0) - m_system->getElectrons().at(j)->getPosition().at(0);
                const double dy = m_system->getElectrons().at(i)->getPosition().at(1) - m_system->getElectrons().at(j)->getPosition().at(1);
                const double dz = m_system->getElectrons().at(i)->getPosition().at(2) - m_system->getElectrons().at(j)->getPosition().at(2);
                const double norm = sqrt(dx*dx + dy*dy + dz*dz);
                electronElectronInteractionEnergy += 1.0 / norm;
            }
        }
        return electronElectronInteractionEnergy;
    } else {
        return 0;
    }
}


Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeLocalEnergy() {
    m_kineticEnergy     = computeKineticEnergy();
    m_potentialEnergy   = computeCoreCorePotentialEnergy() +
                          computeElectronCorePotentialEnergy() +
                          computeElectronElectronPotentialEnergy();
    m_localEnergy = m_kineticEnergy + m_potentialEnergy;
    return m_localEnergy;
}
