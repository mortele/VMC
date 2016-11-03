#include "system.h"
#include "electron.h"
#include "Cores/core.h"
#include "Hamiltonians/hamiltonian.h"

System::System() {
    m_electrons.clear();
    m_cores.clear();
    m_numberOfElectrons = 0;
    m_numberOfDimensions = 3;
}


void System::setHamiltonian(class Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

double System::computeLocalEnergy() {
    return m_hamiltonian->computeLocalEnergy();
}

void System::adjustPositionOfElectron(int       electron,
                                      int       dimension,
                                      double    adjustment) {
    m_electrons.at(electron)->adjustPosition(adjustment, dimension);
}

void System::addCore(Core* core) {
    m_cores.push_back(core);
    for (Electron* electron : core->getElectrons()) {
        m_electrons.push_back(electron);
        m_numberOfElectrons += 1;
    }
}
