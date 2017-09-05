#include "gaussianorbital.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/Orbitals/contractedgaussian.h"


GaussianOrbital::GaussianOrbital(std::string basisFileName) {
    m_parser = new HartreeFockBasisParser();
    m_parser->parseBasisFile(basisFileName);
    m_basis = m_parser->getBasis();
    m_spinUpCoefficients   = m_parser->getSpinUpCoefficients();
    m_spinDownCoefficients = m_parser->getSpinDownCoefficients();
    m_basisSize = m_parser->getBasisSize();
}

double GaussianOrbital::evaluate(double x,
                                 double y,
                                 double z,
                                 int index) {
    double value = 0;
    for (int i=0; i < m_basisSize; i++) {
        value += m_spinUpCoefficients(i,index) * m_basis.at(index)->evaluate(x,y,z);
    }
    return value;
    //return m_basis.at(index)->evaluate(x,y,z);
}

double GaussianOrbital::computeDerivative(double x,
                                          double y,
                                          double z,
                                          int index,
                                          int dimension) {

    if (dimension==0) {
        double value = 0;
        for (int i=0; i < m_basisSize; i++) {
            value += m_spinUpCoefficients(i,index) * m_basis.at(index)->xDerivative(x,y,z);
        }
        return value;
        //return m_basis.at(index)->xDerivative(x,y,z);
    } else if (dimension==1) {
        double value = 0;
        for (int i=0; i < m_basisSize; i++) {
            value += m_spinUpCoefficients(i,index) * m_basis.at(index)->yDerivative(x,y,z);
        }
        return value;
        //return m_basis.at(index)->yDerivative(x,y,z);
    } else {
        double value = 0;
        for (int i=0; i < m_basisSize; i++) {
            value += m_spinUpCoefficients(i,index) * m_basis.at(index)->zDerivative(x,y,z);
        }
        return value;
        //return m_basis.at(index)->zDerivative(x,y,z);
    }
}

double GaussianOrbital::computeLaplacian(double x,
                                         double y,
                                         double z,
                                         int index) {
    double value = 0;
    for (int i=0; i < m_basisSize; i++) {
        value += m_spinUpCoefficients(i,index) * m_basis.at(index)->calculateLaplacian(x,y,z);
    }
    return value;
    //return m_basis.at(index)->calculateLaplacian(x,y,z);
}
