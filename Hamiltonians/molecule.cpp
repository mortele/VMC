#include "molecule.h"
#include "system.h"
#include "Cores/atom.h"
#include <armadillo>


Molecule::Molecule(System* system) :
    Hamiltonian(system) {
}

double Molecule::computeLocalEnergy() {

}

void Molecule::setup() {
    m_atoms.clear();
    m_atoms.reserve(m_system->getCores().size());

    for (Core* atom : m_system->getCores()) {
        m_atoms.push_back(atom);
    }
}
