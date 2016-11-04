#include "Cores/atom.h"
#include "electron.h"
#include "system.h"
#include <iostream>
#include <armadillo>
#include <RandomNumberGenerator/random.h>

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

void Atom::createElectrons() {
    findAtomSize();
    cout << m_size << endl;
    for (int i = 0; i < m_charge; i++) {
        m_electrons.push_back(new Electron(vec {Random::nextGaussian(m_position(0), m_size),
                                                Random::nextGaussian(m_position(1), m_size),
                                                Random::nextGaussian(m_position(2), m_size)},
                                           i % 2));
    }
}

Atom::Atom(System* system, arma::vec position, int charge) :
        Core(system, position) {
    m_charge            = charge;
    m_generalizedCharge = (double) charge;
    m_electrons.clear();
    m_electrons.reserve(charge);
    createElectrons();
}

double Atom::computeCoreCoreInteraction() {
    // Implicit assumption: All the other cores in m_system are atoms (not i.e.
    // quantum dots).

    double interactionEnergy = 0;
    for (Core* core : m_system->getCores()) {
        if (core != this) {
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
    return interactionEnergy;
}

void Atom::findAtomSize() {
    // Info from:
    // https://en.wikipedia.org/wiki/Atomic_radius#Calculated_atomic_radii

    // Find size in pm.
    if        (m_charge == 1) {
        m_size = 53;
    } else if (m_charge == 2) {
        m_size = 31;
    } else if (m_charge == 3) {
        m_size = 167;
    } else if (m_charge == 4) {
        m_size = 112;
    } else if (m_charge == 5) {
        m_size = 87;
    } else if (m_charge == 6) {
        m_size = 67;
    } else if (m_charge == 7) {
        m_size = 56;
    } else if (m_charge == 8) {
        m_size = 48;
    } else if (m_charge == 9) {
        m_size = 42;
    } else if (m_charge == 10) {
        m_size = 38;
    } else if (m_charge == 11) {
        m_size = 190;
    } else if (m_charge == 12) {
        m_size = 145;
    } else if (m_charge == 13) {
        m_size = 118;
    } else if (m_charge == 14) {
        m_size = 111;
    } else if (m_charge == 15) {
        m_size = 98;
    } else if (m_charge == 16) {
        m_size = 88;
    } else if (m_charge == 17) {
        m_size = 79;
    } else if (m_charge == 18) {
        m_size = 71;
    }
    // Convert size to atomic units, Bohr radii (a0).
    const double pmToBohrRadius = 0.0188972613;
    m_size *= pmToBohrRadius;
}















