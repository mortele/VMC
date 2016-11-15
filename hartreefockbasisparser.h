#pragma once
#include <string>
#include <vector>
#include <armadillo>

class HartreeFockBasisParser {
private:
    int             m_numberOfElectrons;
    int             m_numberOfSpinUpElectrons;
    int             m_numberOfSpinDownElectrons;
    int             m_basisSize;
    int             m_numberOfAtoms;
    arma::mat       m_spinUpCoefficients;
    arma::mat       m_spinDownCoefficients;
    std::string     m_fileName;
    std::vector<int>                        m_atomCharges;
    std::vector<arma::vec>                  m_atomPositions;
    std::vector<class ContractedGaussian*>  m_basis;

public:
    void parseBasisFile(std::string fileName);

    int getNumberOfElectrons()               { return m_numberOfElectrons; }
    int getNumberOfSpinUpElectrons()         { return m_numberOfSpinUpElectrons; }
    int getNumberOfSpinDownElectrons()       { return m_numberOfSpinDownElectrons; }
    int getBasisSize()                       { return m_basisSize; }
    int getNumberOfAtoms()                   { return m_numberOfAtoms; }
    arma::mat getSpinUpCoefficients()        { return m_spinUpCoefficients; }
    arma::mat getSpinDownCoefficients()      { return m_spinDownCoefficients; }
    std::vector<int>       getAtomCharges()  { return m_atomCharges; }
    std::vector<arma::vec> getAtomPosition() { return m_atomPositions; }
    std::vector<class ContractedGaussian*> getBasis() { return m_basis; }
};

