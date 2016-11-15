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

PrimitiveGaussian::PrimitiveGaussian(int    x,
                                     int    y,
                                     int    z,
                                     double alpha,
                                     double constant) :
        m_x         (x),
        m_y         (x),
        m_z         (z),
        m_alpha     (alpha),
        m_constant  (constant) {
}

double PrimitiveGaussian::operator()(double x, double y, double z) {
    const double r2 = x*x + y*y + z*z;
    return m_constant * pow(x, m_x) * pow(y, m_y) * pow(z, m_z) * exp(- m_alpha * r2);
}

double PrimitiveGaussian::xDerivative(double x, double y, double z) {

}

