#include "hamiltonian.h"
#include "WaveFunctions/hydrogenwavefunction.h"
#include "Cores/core.h"
#include "electron.h"
#include "system.h"
#include <armadillo>


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
        for (int i = 0; i < m_system->getElectrons().size(); i++) {
            for (int j = i+1; j < m_system->getElectrons().size(); j++) {
                vec dr = m_system->getElectrons().at(i)->getPosition() -
                         m_system->getElectrons().at(j)->getPosition();
                electronElectronInteractionEnergy += 1.0 / norm(dr);
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

    double kinetic = computeKineticEnergy();
    double corecore = computeCoreCorePotentialEnergy();
    double coreelectron = computeElectronCorePotentialEnergy();
    double electronelectron = computeElectronElectronPotentialEnergy();
    double tot = kinetic + corecore + coreelectron + electronelectron;

    const vec position1 = m_system->getElectrons().at(0)->getPosition();
    const vec position2 = m_system->getElectrons().at(1)->getPosition();
    const double r1  = norm(position1);
    const double r2  = norm(position2);
    const double r12        = norm(position1 - position2);
    const double r12Inverse = 1.0 / r12;
    const double r1Dotr2    = dot(position1, position2);
    const double r1Inverse  = 1.0 / r1;
    const double r2Inverse  = 1.0 / r2;
    const double Z = 2.0;
    double a = 1.843;
    double b = 0.347;
    double b1 = b * r12;
    double b2 = 1 + b1;
    double b3 = 1/(2 * b2 * b2);
    double prikk = r1Dotr2;
    double E_L1 = (a - Z) * (1 / r1  + 1 / r2) + 1 / r12 - a*a;
    double E_L2 = E_L1  + b3 * ((a * (r1 + r2)) / (r12)  * (1  - (prikk / (r1 * r2))) - b3 - 2 / r12 + ((2*b) / b2)); //

    //return E_L2;    //std::cout << E_L2 - kinetic << std::endl;

    m_localEnergy = computeKineticEnergy() +
                    computeCoreCorePotentialEnergy() +
                    computeElectronCorePotentialEnergy() +
                    computeElectronElectronPotentialEnergy();
    return m_localEnergy;
}
