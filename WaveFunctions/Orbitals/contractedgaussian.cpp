#include "contractedgaussian.h"
#include "WaveFunctions/Orbitals/primitivegaussian.h"


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

double ContractedGaussian::xDerivative(double x, double y, double z) {
    double result = 0;
    for (int i = 0; i < m_numberOfPrimitives; i++) {
        result += m_primitives.at(i)->xDerivative(x,y,z);
    }
    return result;
}

double ContractedGaussian::yDerivative(double x, double y, double z) {
    double result = 0;
    for (int i = 0; i < m_numberOfPrimitives; i++) {
        result += m_primitives.at(i)->yDerivative(x,y,z);
    }
    return result;
}

double ContractedGaussian::zDerivative(double x, double y, double z) {
    double result = 0;
    for (int i = 0; i < m_numberOfPrimitives; i++) {
        result += m_primitives.at(i)->zDerivative(x,y,z);
    }
    return result;
}

double ContractedGaussian::calculateLaplacian(double x, double y, double z) {
    const double xA         = x - m_x;
    const double yA         = y - m_y;
    const double zA         = z - m_z;
    double totalValue       = 0;
    double totalLaplacian   = 0;
    for (PrimitiveGaussian* primitive : m_primitives) {
        const double value   = primitive->getCurrentValue(); //(*primitive)(xA,yA,zA);
        totalValue          += value;
        const double ddx     = primitive->xxDerivative(xA,yA,zA);
        const double ddy     = primitive->yyDerivative(xA,yA,zA);
        const double ddz     = primitive->zzDerivative(xA,yA,zA);
        totalLaplacian      += (ddx + ddy + ddz) * value;
    }
    m_currentValue = totalValue;
    return totalLaplacian;
}








