#pragma once

class Sampler {
private:
    class System*       m_system;
    class Hamiltonian*  m_hamiltonian;



public:
    Sampler(class System* system);
    void setup();
};
