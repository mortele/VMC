#pragma once
#include <vector>
#include <armadillo>


class Core {
protected:
    double                          m_generalizedCharge;
    double                          m_size;
    arma::vec                       m_position;
    class System*                   m_system;
    std::vector<class Electron*>    m_electrons;

public:
    Core(class System* system, arma::vec position);
    virtual double computeCoreCoreInteraction() { return 0; }
    virtual double computeElectronCoreInteraction() = 0;
    virtual std::string getInfo() { return "Unknown"; }

    double      getGeneralizedCharge         () { return m_generalizedCharge; }
    double      getSize                      () { return m_size; }
    arma::vec   getPosition                  () { return m_position; }
    std::vector<class Electron*> getElectrons() { return m_electrons; }
};

