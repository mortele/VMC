#include "primitivegaussian.h"
#include <cmath>

using std::exp;

double PrimitiveGaussian::pow(double a, int n) {
    double result = 1;
    for (int i=0; i<n; i++) {
        result *= a;
    }
    return result;
}

PrimitiveGaussian::PrimitiveGaussian(int    i,
                                     int    j,
                                     int    k,
                                     double Ax,
                                     double Ay,
                                     double Az,
                                     double alpha,
                                     double constant) :
        m_x         (i),
        m_y         (j),
        m_z         (k),
        m_Ax        (Ax),
        m_Ay        (Ay),
        m_Az        (Az),
        m_alpha     (alpha),
        m_constant  (constant) {
}

double PrimitiveGaussian::operator()(double x, double y, double z) {
    const double r2 = x*x + y*y + z*z;
    const double value = m_constant * pow(x, m_x) * pow(y, m_y) * pow(z, m_z) * exp(- m_alpha * r2);
    m_currentValue = value;
    return value;
}

double PrimitiveGaussian::xDerivative(double x, double y, double z) {
    //const double value = (*this)(x, y, z);
    return m_currentValue * ((m_x == 0) ? - 2 * m_alpha * x : m_x / x - 2 * m_alpha * x);
}

double PrimitiveGaussian::yDerivative(double x, double y, double z) {
    //const double value = (*this)(x, y, z);
    return m_currentValue * ((m_y == 0) ? - 2 * m_alpha * y : m_y / y - 2 * m_alpha * y);
}

double PrimitiveGaussian::zDerivative(double x, double y, double z) {
    //const double value = (*this)(x, y, z);
    return m_currentValue * ((m_z == 0) ? - 2 * m_alpha * z : m_z / z - 2 * m_alpha * z);
}

double PrimitiveGaussian::xxDerivative(double x, double y, double z) {
    //m_currentValue = (*this)(x,y,z);
    const double a2 = 2*m_alpha;
    if (m_x==0) {
        return a2*(a2*x*x - 1);
    } else if (m_x==1) {
        return a2*(a2*x*x - 3);
    } else {
        return (a2*x*x*(a2*x*x - 2*m_x - 1) + m_x*m_x - m_x) / (x*x);
    }
}

double PrimitiveGaussian::yyDerivative(double x, double y, double z) {
    //const double a2 = 2*m_alpha;
    if (m_y==0) {
        return a2*(a2*y*y - 1);
    } else if (m_y==1) {
        return a2*(a2*y*y - 3);
    } else {
        return (a2*y*y*(a2*y*y - 2*m_y - 1) + m_y*m_y - m_y) / (y*y);
    }
}

double PrimitiveGaussian::zzDerivative(double x, double y, double z) {
    //const double a2 = 2*m_alpha;
    if (m_z==0) {
        return a2*(a2*z*z - 1);
    } else if (m_z==1) {
        return a2*(a2*z*z - 3);
    } else {
        return (a2*z*z*(a2*z*z - 2*m_z - 1) + m_z*m_z - m_z) / (z*z);
    }
}

