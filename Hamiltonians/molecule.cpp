#include "molecule.h"
#include "system.h"
#include "atom.h"
#include <armadillo>


Molecule::Molecule(System* system) :
    Hamiltonian(system) {
}

double Molecule::computeLocalEnergy() {

}

void Molecule::setup() {
    m_atoms.clear();
    m_atoms.reserve(m_system->getAtoms().size());

    for (Atom* atom : m_system->getAtoms()) {
        m_atoms.push_back(atom);
    }
}
