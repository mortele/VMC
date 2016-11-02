#pragma once
#include <armadillo>

class Electron {
private:
    int       m_spin;
    arma::vec m_position;

public:
    Electron(arma::vec position, int spin);

    void setPosition    (arma::vec position);
    void setPosition    (double position, int dimension);
    void adjustPosition (double adjustment, int dimension);

    int         getSpin()       { return m_spin; }
    arma::vec   getPosition()   { return m_position; }
};

