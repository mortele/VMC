#pragma once
#include "Hamiltonians/hamiltonian.h"
#include <vector>


class Molecule : public Hamiltonian {
private:
    class System*               m_system;
    std::vector<class Core*>    m_atoms;

public:
    Molecule(class System* system);
    double computeLocalEnergy();
    void setup();
};

