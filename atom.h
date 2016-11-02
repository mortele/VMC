#pragma once
#include <vector>
#include <armadillo>

class Atom {
private:
    int                             m_charge;
    arma::vec                       m_position;
    std::vector<class Electron*>    m_electrons;

public:
    Atom(arma::vec position, double charge);

    std::vector<class Electron*> getElectrons() { return m_electrons; }
};

