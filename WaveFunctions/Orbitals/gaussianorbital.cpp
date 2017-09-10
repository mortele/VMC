#include "gaussianorbital.h"
#include "hartreefockbasisparser.h"
#include "WaveFunctions/Orbitals/contractedgaussian.h"
#include <iostream>
#include <iomanip>
#include "system.h"

using std::cout;
using std::endl;
using std::setprecision;

GaussianOrbital::GaussianOrbital(System* system, std::string basisFileName) {
    m_parser = new HartreeFockBasisParser();
    m_parser->parseBasisFile(system, basisFileName);
    m_basis = m_parser->getBasis();
    m_spinUpCoefficients   = m_parser->getSpinUpCoefficients();
    m_spinDownCoefficients = m_parser->getSpinDownCoefficients();
    m_basisSize = m_parser->getBasisSize();
}

double GaussianOrbital::evaluate(double x,
                                 double y,
                                 double z,
                                 int index,
                                 int spin) {
    double value = 0;
    cout << "x,y,z=" << x << " " << y << " " << z << endl;
    for (int i=0; i < m_basisSize; i++) {
        double c = (spin==1 ? m_spinUpCoefficients(i,index) : m_spinDownCoefficients(i,index));
        if (fabs(c) > 1e-10) {
            cout << "index,spin,i: " << index << "," << spin << "," << i << " " << c << " " << m_basis.at(i)->evaluate(x,y,z) << " " <<  c*m_basis.at(i)->evaluate(x,y,z) << endl;
            cout << *m_basis.at(i) << endl;
        }
        value += (spin==1 ? m_spinUpCoefficients(i,index) : m_spinDownCoefficients(i,index))
                 * m_basis.at(i)->evaluate(x,y,z);
    }
    cout << "TOTAL: " << value << endl << endl;
    return value;
}

double GaussianOrbital::computeDerivative(double x,
                                          double y,
                                          double z,
                                          int index,
                                          int dimension,
                                          int spin) {
    if (dimension==0) {
        double value = 0;
        for (int i=0; i < m_basisSize; i++) {
            value += (spin==1 ? m_spinUpCoefficients(i,index) : m_spinDownCoefficients(i,index))
                    * m_basis.at(i)->xDerivative(x,y,z);
        }
        return value;
    } else if (dimension==1) {
        double value = 0;
        for (int i=0; i < m_basisSize; i++) {
            value += (spin==1 ? m_spinUpCoefficients(i,index) : m_spinDownCoefficients(i,index))
                    * m_basis.at(i)->yDerivative(x,y,z);
        }
        return value;
    } else {
        double value = 0;
        for (int i=0; i < m_basisSize; i++) {
            value += (spin==1 ? m_spinUpCoefficients(i,index) : m_spinDownCoefficients(i,index))
                    * m_basis.at(i)->zDerivative(x,y,z);
        }
        return value;
    }
}

double GaussianOrbital::computeLaplacian(double x,
                                         double y,
                                         double z,
                                         int index,
                                         int spin) {
    double value = 0;
    for (int i=0; i < m_basisSize; i++) {
        value += (spin==1 ? m_spinUpCoefficients(i,index) : m_spinDownCoefficients(i,index))
                * m_basis.at(i)->calculateLaplacian(x,y,z);
    }
    return value;
}
