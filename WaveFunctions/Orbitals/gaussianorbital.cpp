#include "gaussianorbital.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/Orbitals/contractedgaussian.h"


GaussianOrbital::GaussianOrbital(std::string basisFileName) {
    m_parser = new HartreeFockBasisParser();
    m_parser->parseBasisFile(basisFileName);
    m_basis = m_parser->getBasis();
}

double GaussianOrbital::evaluate(double x,
                                 double y,
                                 double z,
                                 int index) {
    return m_basis.at(index)->evaluate(x,y,z);
}

double GaussianOrbital::computeDerivative(double x,
                                          double y,
                                          double z,
                                          int index,
                                          int dimension) {
    if (dimension==0) {
        return m_basis.at(index)->xDerivative(x,y,z);
    } else if (dimension==1) {
        return m_basis.at(index)->yDerivative(x,y,z);
    } else {
        return m_basis.at(index)->zDerivative(x,y,z);
    }
}

double GaussianOrbital::computeLaplacian(double x,
                                         double y,
                                         double z,
                                         int index) {
    return m_basis.at(index)->calculateLaplacian(x,y,z);
}
