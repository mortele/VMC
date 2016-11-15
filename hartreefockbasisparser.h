#pragma once
#include <string>
#include <vector>
#include <armadillo>

class HartreeFockBasisParser {
private:
    int         m_numberOfElectrons;
    int         m_numberOfSpinUpElectrons;
    int         m_numberOfSpinDownElectrons;
    int         m_basisSize;
    int         m_numberOfAtoms;
    std::string m_fileName;
    std::vector<arma::vec> m_atomPositions;

public:
    class WaveFunction* parseBasisFile(std::string fileName);
};

