#pragma once
#include "Cores/core.h"
#include <vector>
#include <armadillo>

class Atom : public Core {
private:
    int m_charge;

public:
    Atom(class System* system, arma::vec position, int charge);

    double computeCoreCoreInteraction();
    double computeElectronCoreInteraction();
};

