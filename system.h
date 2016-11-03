#pragma once
#include <vector>

class System {
private:
    int                             m_numberOfDimensions = 3;
    int                             m_numberOfElectrons;
    class Hamiltonian*              m_hamiltonian;
    std::vector<class Core*>        m_cores;
    std::vector<class Electron*>    m_electrons;

public:
    System();

    double computeLocalEnergy();

    void adjustPositionOfElectron(int electron, int dimension, double adjustment);
    void addCore        (class Core* core);
    void setHamiltonian (class Hamiltonian* hamiltonian);

    int getNumberOfDimensions () { return m_numberOfDimensions; }
    int getNumberOfElectrons  () { return m_numberOfElectrons; }
    class Hamiltonian*           getHamiltonian() { return m_hamiltonian; }
    std::vector<class Core*>     getCores()       { return m_cores; }
    std::vector<class Electron*> getElectrons()   { return m_electrons; }
};

