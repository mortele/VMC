#pragma once

class Hamiltonian {
private:
    class System* m_system;

public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy() = 0;
};

