#include "contractedgaussian.h"
#include "WaveFunctions/Orbitals/primitivegaussian.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;
using std::setprecision;


ContractedGaussian::ContractedGaussian(int      numberOfPrimitives,
                                       double   x,
                                       double   y,
                                       double   z) {
    m_numberOfPrimitives = numberOfPrimitives;
    m_primitives.reserve(numberOfPrimitives);
    m_x = x;
    m_y = y;
    m_z = z;
}

void ContractedGaussian::addPrimitive(PrimitiveGaussian* primitive) {
    m_primitives.push_back(primitive);
}

double ContractedGaussian::operator()(double x, double y, double z) {
    double result = 0;
    for (PrimitiveGaussian* primitive : m_primitives) {
        result += (*primitive)(x - m_x, y - m_y, z - m_z);
    }
    m_currentValue = result;
    return result;
}

double ContractedGaussian::evaluate(double x, double y, double z) {
    return (*this)(x,y,z);
}

double ContractedGaussian::xDerivative(double x, double y, double z) {
    evaluate(x,y,z);
    double result = 0;
    for (int i = 0; i < m_numberOfPrimitives; i++) {
        PrimitiveGaussian*& primitive = m_primitives.at(i);
        result += primitive->xDerivative(x,y,z);
    }
    return result;
}

double ContractedGaussian::yDerivative(double x, double y, double z) {
    double result = 0;
    for (int i = 0; i < m_numberOfPrimitives; i++) {
        PrimitiveGaussian*& primitive = m_primitives.at(i);
        result += primitive->yDerivative(x,y,z);
    }
    return result;
}

double ContractedGaussian::zDerivative(double x, double y, double z) {
    double result = 0;
    for (int i = 0; i < m_numberOfPrimitives; i++) {
        PrimitiveGaussian*& primitive = m_primitives.at(i);
        result += primitive->zDerivative(x,y,z);
    }
    return result;
}

double ContractedGaussian::calculateLaplacian(double x, double y, double z) {
    evaluate(x,y,z);
    const double xA         = x - m_x;
    const double yA         = y - m_y;
    const double zA         = z - m_z;
    double totalLaplacian   = 0;
    for (PrimitiveGaussian* primitive : m_primitives) {
        const double ddx     = primitive->xxDerivative(xA,yA,zA);
        const double ddy     = primitive->yyDerivative(xA,yA,zA);
        const double ddz     = primitive->zzDerivative(xA,yA,zA);
        totalLaplacian      += (ddx + ddy + ddz) * primitive->getCurrentValue();
    }
    return totalLaplacian;
}

std::ostream& operator<<(std::ostream& stream, const ContractedGaussian& contracted) {

    cout << "nucleus: " << "(" << contracted.m_x << ","
                               << contracted.m_y << ","
                               << contracted.m_z << ")" << endl;

    int i = 0;
    for (PrimitiveGaussian* primitive : contracted.m_primitives) {
        stream << "   * primitive" << i << ": " << *primitive << endl;
        i++;
    }
    return stream;
}







