#pragma once
#include "Cores/core.h"
#include <armadillo>

class HarmonicOscillator : public Core {
private:
    int m_numberOfElectrons;

    void createElectrons();

public:
    HarmonicOscillator(class System*    system,
                       arma::vec        position,
                       int              numberOfElectrons,
                       double           omega);
    double computeElectronCoreInteraction();
};
