#include "primitivegaussian.h"
#include <cmath>
#include <iomanip>

using std::exp;

inline double PrimitiveGaussian::pow(double a, int n) {
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
        m_i         (i),
        m_j         (j),
        m_k         (k),
        m_Ax        (Ax),
        m_Ay        (Ay),
        m_Az        (Az),
        m_alpha     (alpha),
        m_constant  (constant) {
}

double PrimitiveGaussian::operator()(double x, double y, double z) {
    const double r2 = x*x + y*y + z*z;
    const double value = m_constant * pow(x, m_i) * pow(y, m_j) * pow(z, m_k) * exp(- m_alpha * r2);
    m_currentValue = value;
    return value;
}

double PrimitiveGaussian::xDerivative(double x, double , double ) {
    double value = m_currentValue;
    return value * ((m_i == 0) ? - 2 * m_alpha * x : m_i / x - 2 * m_alpha * x);
}

double PrimitiveGaussian::yDerivative(double , double y, double ) {
    double value = m_currentValue;
    return value * ((m_j == 0) ? - 2 * m_alpha * y : m_j / y - 2 * m_alpha * y);
}

double PrimitiveGaussian::zDerivative(double , double , double z) {
    double value = m_currentValue;
    return value * ((m_k == 0) ? - 2 * m_alpha * z : m_k / z - 2 * m_alpha * z);
}

double PrimitiveGaussian::xxDerivative(double x, double , double ) {
    const double a2 = 2*m_alpha;
    if (m_i==0) {
        return a2*(a2*x*x - 1);
    } else if (m_i==1) {
        return a2*(a2*x*x - 3);
    } else {
        return (a2*x*x*(a2*x*x - 2*m_i - 1) + m_i*m_i - m_i) / (x*x);
    }
}

double PrimitiveGaussian::yyDerivative(double , double y, double ) {
    const double a2 = 2*m_alpha;
    if (m_j==0) {
        return a2*(a2*y*y - 1);
    } else if (m_j==1) {
        return a2*(a2*y*y - 3);
    } else {
        return (a2*y*y*(a2*y*y - 2*m_j - 1) + m_j*m_j - m_j) / (y*y);
    }
}

double PrimitiveGaussian::zzDerivative(double , double , double z) {
    const double a2 = 2*m_alpha;
    if (m_k==0) {
        return a2*(a2*z*z - 1);
    } else if (m_k==1) {
        return a2*(a2*z*z - 3);
    } else {
        return (a2*z*z*(a2*z*z - 2*m_k - 1) + m_k*m_k - m_k) / (z*z);
    }
}


std::ostream& operator<<(std::ostream& stream, const PrimitiveGaussian& primitive) {
    stream << std::setprecision(5) << primitive.m_constant;
    stream << "(";
    stream << primitive.m_i << ",";
    stream << primitive.m_j << ",";
    stream << primitive.m_k << ")";
    stream << " exp(- " << primitive.m_alpha << " r^2)";
    return stream;
}
