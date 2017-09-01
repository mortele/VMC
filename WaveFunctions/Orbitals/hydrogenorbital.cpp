#include "hydrogenorbital.h"

using std::exp;


HydrogenOrbital::HydrogenOrbital(double alpha) {
    m_alpha  = alpha;
    m_alpha2 = m_alpha*m_alpha;
}

double HydrogenOrbital::evaluate1s(double r) {
    return exp(-m_alpha * r);
}

double HydrogenOrbital::evaluate2s(double r) {
    return (1 - m_alpha * 0.5*distance) * exp(-m_alpha * 0.5 * distance);
}

double HydrogenOrbital::evaluate2p(double r, double x) {
    // x ≡ x or y or z, for 2px, 2py, and 2pz, respectively.
    return x*exp(-m_alpha * 0.5 * r);
}

double HydrogenOrbital::computeDerivative1s(double r, double x) {
    // x ≡ x or y or z, the variable we are differentiation w.r.t.
    return -m_alpha * x * exp(-m_alpha*r) / r;
}

double HydrogenOrbital::computeDerivative2s(double r, double x) {
    return 0.25*m_alpha * (m_alpha * r - 4) * x * exp(-0.5* m_alpha*r) / r;
}

double HydrogenOrbital::computeDerivative2px(double  r,
                                             double  x,
                                             double  y,
                                             double  z,
                                             int     i) {
    // i determines which variable we are differentiation w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return -(m_alpha*x*x - 2*r) * exp(-0.5*m_alpha*r) / (2*r);
    } else if (i==1) { // y derivative
        return -m_alpha*x*y * exp(-0.5*m_alpha*r) / (2*r);
    } else { // z derivative
        return -m_alpha*x*z * exp(-0.5*m_alpha*r) / (2*r);
    }
}

double HydrogenOrbital::computeDerivative2py(double r,
                                             double x,
                                             double y,
                                             double z,
                                             int    i) {
    // i determines which variable we are differentiation w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return -m_alpha*x*y * exp(-0.5*m_alpha*r) / (2*r);
    } else if (i==1) { // y derivative
        return -(m_alpha*y*y - 2*r) * exp(-0.5*m_alpha*r) / (2*r);
    } else { // z derivative
        return -m_alpha*y*z * exp(-0.5*m_alpha*r) / (2*r);
    }
}

double HydrogenOrbital::computeDerivative2pz(double r,
                                             double x,
                                             double y,
                                             double z,
                                             int i) {
    // i determines which variable we are differentiation w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return -m_alpha*x*z * exp(-0.5*m_alpha*r) / (2*r);
    } else if (i==1) { // y derivative
        return -m_alpha*y*z * exp(-0.5*m_alpha*r) / (2*r);
    } else { // z derivative
        return -(m_alpha*z*z - 2*r) * exp(-0.5*m_alpha*r) / (2*r);
    }
}

double HydrogenOrbital::computeLaplacian1s(double r) {
    return (m_alpha2 - 2 * m_alpha / r) * exp(-m_alpha * r);
}

double HydrogenOrbital::computeLaplacian2s(double r) {
    return (1.25 * m_alpha2 - 2 * m_alpha / distance - 0.125 * m_alpha2*m_alpha * r) * exp(-0.5* m_alpha * r);
}

double HydrogenOrbital::computeLaplacian2p(double r, double x) {
    // x ≡ x or y or z, for 2px, 2py, and 2pz, respectively.
    return m_alpha*x * (m_alpha*r - 8.) * exp(-0.5*m_alpha*r) / (4*r);
}



























