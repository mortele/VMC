#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/hydrogenwavefunction.h"
#include "Cores/core.h"
#include "electron.h"
#include "system.h"
#include <armadillo>

using arma::vec;
using arma::norm;


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
    double electronElectronInteractionEnergy = 0;
    for (Electron* electron1 : m_system->getElectrons()) {
        for (Electron* electron2 : m_system->getElectrons()) {
            if (electron1 != electron2) { // TODO: Is this horribly bad practice?
                vec     dr          = electron1->getPosition() -
                                      electron2->getPosition();
                electronElectronInteractionEnergy += 1.0 / norm(dr);
            }
        }
    }
    return 0.5 * electronElectronInteractionEnergy;
}


Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeLocalEnergy() {

}
