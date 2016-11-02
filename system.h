#pragma once
#include <vector>

class System {
private:
    class Hamiltonian*              m_hamiltonian;
    std::vector<class Atom*>        m_atoms;
    std::vector<class Electron*>    m_electrons;

public:
    System();
    void setHamiltonian(class Hamiltonian* hamiltonian);
    void addAtom(class Atom* atom);
};

