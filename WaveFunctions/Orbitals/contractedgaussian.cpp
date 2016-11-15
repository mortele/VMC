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
    return result;
}
