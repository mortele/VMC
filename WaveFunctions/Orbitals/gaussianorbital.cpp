#include "gaussianorbital.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/Orbitals/contractedgaussian.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setprecision;

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
        value += m_spinUpCoefficients(i,index) * m_basis.at(i)->evaluate(x,y,z);
        //cout << i << " " << m_spinUpCoefficients(i,index) * m_basis.at(i)->evaluate(x,y,z) << " " << value <<  endl;
    }
    return value;
    //return m_basis.at(index)->evaluate(x,y,z);
}

double GaussianOrbital::computeDerivative(double x,
                                          double y,
                                          double z,
                                          int index,
                                          int dimension) {
    bool debugprint = false;
    if (debugprint) cout << x << " " << y << " " << z << endl;
    if (debugprint) cout << "=======================" << endl;
    if (debugprint) cout << std::sqrt(x*x+y*y+z*z) << endl;
    if (debugprint) cout << "=======================" << endl;

    if (dimension==0) {
        double value = 0;
        for (int i=0; i < m_basisSize; i++) {
            if (debugprint) cout << "C " << m_spinUpCoefficients(i,index) << endl;
            value += m_spinUpCoefficients(i,index) * m_basis.at(i)->xDerivative(x,y,z);
        }
        if (debugprint) cout << "total val: " << value << endl;
        if (debugprint) exit(1);

        return value;
        //return m_basis.at(index)->xDerivative(x,y,z);
    } else if (dimension==1) {
        double value = 0;
        for (int i=0; i < m_basisSize; i++) {
            value += m_spinUpCoefficients(i,index) * m_basis.at(i)->yDerivative(x,y,z);
        }
        return value;
        //return m_basis.at(index)->yDerivative(x,y,z);
    } else {
        double value = 0;
        for (int i=0; i < m_basisSize; i++) {
            value += m_spinUpCoefficients(i,index) * m_basis.at(i)->zDerivative(x,y,z);
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
        value += m_spinUpCoefficients(i,index) * m_basis.at(i)->calculateLaplacian(x,y,z);
    }
    return value;
    //return m_basis.at(index)->calculateLaplacian(x,y,z);
}
