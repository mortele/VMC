#include "Cores/atom.h"
#include "electron.h"
#include "system.h"
#include <iostream>

using std::cout;
using std::endl;

Atom::Atom(System* system, arma::vec position, int charge) :
        Core(system, position) {
    m_charge            = charge;
    m_generalizedCharge = (double) charge;
    m_electrons.clear();
    m_electrons.reserve(charge);
}

double Atom::computeCoreCoreInteraction() {
    // Implicit assumption: All the other cores in m_system are atoms (not i.e.
    // quantum dots).

    double interactionEnergy = 0;
    for (Core* core : m_system->getCores()) {
        if (core != this) { // TODO : Is this horribly bad practice?
            arma::vec   dr              = core->getPosition() - m_position;
            double      rInverse        = 1.0 / norm(dr);
            double      chargeProduct   = m_charge * core->getGeneralizedCharge();
            interactionEnergy          += chargeProduct * rInverse;
        }
    }
    return interactionEnergy;
}

double Atom::computeElectronCoreInteraction() {
    double interactionEnergy = 0;
    for (Electron* electron : m_system->getElectrons()) {
        arma::vec   dr          = electron->getPosition() - m_position;
        double      rInverse    = 1.0 / norm(dr);
        interactionEnergy      -= m_charge * rInverse;
    }

    return 0;
}
