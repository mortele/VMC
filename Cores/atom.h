#pragma once
#include "Cores/core.h"
#include <vector>
#include <string>
#include <armadillo>

class Atom : public Core {
private:
    int     m_charge;
    std::string m_atomName;

    void createElectrons();
    void findAtomSize();
    std::string findAtomName();

public:
    Atom(class System* system, arma::vec position, int charge);
    double computeCoreCoreInteraction();
    double computeElectronCoreInteraction();
    std::string getInfo();
};

