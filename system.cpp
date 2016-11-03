#include "system.h"
#include "electron.h"
#include "atom.h"
#include "Hamiltonians/hamiltonian.h"

System::System() {
    m_electrons.clear();
    m_atoms.clear();
    m_numberOfElectrons = 0;
    m_numberOfDimensions = 3;
}


void System::setHamiltonian(class Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::addAtom(Atom* atom) {
    m_atoms.push_back(atom);
    for (Electron* electron : atom->getElectrons()) {
        m_electrons.push_back(electron);
        m_numberOfElectrons += 1;
    }
}

double System::computeLocalEnergy() {
    return m_hamiltonian->computeLocalEnergy();
}

void System::adjustPositionOfElectron(int       electron,
                                      int       dimension,
                                      double    adjustment) {
    m_electrons.at(electron)->adjustPosition(adjustment, dimension);
}
