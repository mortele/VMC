#include "system.h"
#include "electron.h"
#include "atom.h"
#include "Hamiltonians/hamiltonian.h"

System::System() {
    m_electrons.clear();
    m_atoms.clear();
}


void System::setHamiltonian(class Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::addAtom(Atom* atom) {
    m_atoms.push_back(atom);
    for (Electron* electron : atom->getElectrons()) {
        m_electrons.push_back(electron);
    }
}
