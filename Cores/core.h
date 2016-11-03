#pragma once
#include <vector>
#include <armadillo>


class Core {
protected:
    arma::vec                       m_position;
    std::vector<class Electron*>    m_electrons;

public:
    Core(arma::vec position);
};

