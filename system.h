#pragma once
#include <vector>

class System {
private:
    int                             m_numberOfDimensions = 3;
    int                             m_numberOfElectrons;
    class Hamiltonian*              m_hamiltonian;
    std::vector<class Atom*>        m_atoms;
    std::vector<class Electron*>    m_electrons;

public:
    System();

    double computeLocalEnergy();

    void adjustPositionOfElectron(int electron, int dimension, double adjustment);
    void addAtom        (class Atom* atom);
    void setHamiltonian (class Hamiltonian* hamiltonian);

    int getNumberOfDimensions () { return m_numberOfDimensions; }
    int getNumberOfElectrons  () { return m_numberOfElectrons; }
    class Hamiltonian*           getHamiltonian() { return m_hamiltonian; }
    std::vector<class Atom*>     getAtoms()       { return m_atoms; }
    std::vector<class Electron*> getElectrons()   { return m_electrons; }
};

