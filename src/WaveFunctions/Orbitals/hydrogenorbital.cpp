#include "hydrogenorbital.h"
#include <iostream>
#include <cmath>

using std::exp;
using std::cout;
using std::endl;
using std::acos;
using std::sqrt;

HydrogenOrbital::HydrogenOrbital(double alpha) :
        Orbital() {

    double pi = acos(-1);

    m_alpha  = alpha;
    m_alpha2 = m_alpha*m_alpha;
    double a3 = pow(m_alpha,3);
    double a5 = pow(m_alpha,5);
    m_1sNormalization = sqrt(a3/pi);
    m_2sNormalization = sqrt(a3/(8*pi));
    m_2pNormalization = sqrt(a5/(32*pi));
}


double HydrogenOrbital::evaluate(double x, double y, double z, int index, int) {
    const double r = sqrt(x*x + y*y + z*z);
    if (index==0) {
        return evaluate1s(r);
    } else if (index==1) {
        return evaluate2s(r);
    } else if (index==2) {
        return evaluate2p(r, x);
    } else if (index==3) {
        return evaluate2p(r, y);
    } else if (index==4) {
        return evaluate2p(r, z);
    } else {
        cout << "0 Invalid hydrogen orbital index: " << index << endl;
        exit(1);
    }
}

double HydrogenOrbital::computeDerivativeX(double x, double y, double z, int index) {
    const double r = sqrt(x*x + y*y + z*z);
    if (index==0) {
        return computeDerivative1s(r, x);
    } else if (index==1) {
        return computeDerivative2s(r, x);
    } else if (index==2) {
        return computeDerivative2px(r, x, y, z, 0);
    } else if (index==3) {
        return computeDerivative2py(r, x, y, z, 0);
    } else if (index==4) {
        return computeDerivative2pz(r, x, y, z, 0);
    } else {
        cout << "1 Invalid hydrogen orbital index: " << index << endl;
        exit(1);
    }
}

double HydrogenOrbital::computeDerivativeY(double x, double y, double z, int index) {
    const double r = sqrt(x*x + y*y + z*z);
    if (index==0) {
        return computeDerivative1s(r, y);
    } else if (index==1) {
        return computeDerivative2s(r, y);
    } else if (index==2) {
        return computeDerivative2px(r, x, y, z, 1);
    } else if (index==3) {
        return computeDerivative2py(r, x, y, z, 1);
    } else if (index==4) {
        return computeDerivative2pz(r, x, y, z, 1);
    } else {
        cout << "2 Invalid hydrogen orbital index: " << index << endl;
        exit(1);
    }
}

double HydrogenOrbital::computeDerivativeZ(double x, double y, double z, int index) {
    const double r = sqrt(x*x + y*y + z*z);
    if (index==0) {
        return computeDerivative1s(r, z);
    } else if (index==1) {
        return computeDerivative2s(r, z);
    } else if (index==2) {
        return computeDerivative2px(r, x, y, z, 2);
    } else if (index==3) {
        return computeDerivative2py(r, x, y, z, 2);
    } else if (index==4) {
        return computeDerivative2pz(r, x, y, z, 2);
    } else {
        cout << "3 Invalid hydrogen orbital index: " << index << endl;
        exit(1);
    }
}

double HydrogenOrbital::computeDerivative(double x,
                                          double y,
                                          double z,
                                          int    index,
                                          int    dimension,
                                          int    ) {
    if (dimension==0) {
        return computeDerivativeX(x,y,z,index);
    } else if (dimension==1) {
        return computeDerivativeY(x,y,z,index);
    } else {
        return computeDerivativeZ(x,y,z,index);
    }
}

double HydrogenOrbital::computeLaplacian(double x, double y, double z, int index, int) {
    const double r = sqrt(x*x + y*y + z*z);
    if (index==0) {
        return computeLaplacian1s(r);
    } else if (index==1) {
        return computeLaplacian2s(r);
    } else if (index==2) {
        return computeLaplacian2p(r, x);
    } else if (index==3) {
        return computeLaplacian2p(r, y);
    } else if (index==4) {
        return computeLaplacian2p(r, z);
    } else {
        cout << "4 Invalid hydrogen orbital index: " << index << endl;
        exit(1);
    }
}

double HydrogenOrbital::evaluate1s(double r) {
    return m_1sNormalization * exp(-m_alpha * r);
}

double HydrogenOrbital::evaluate2s(double r) {
    return m_2sNormalization * (1 - m_alpha * 0.5*r) * exp(-m_alpha * 0.5 * r);
}

double HydrogenOrbital::evaluate2p(double r, double x) {
    // x ≡ x or y or z, for 2px, 2py, and 2pz, respectively.
    return m_2pNormalization * x*exp(-m_alpha * 0.5 * r);
}

double HydrogenOrbital::computeDerivative1s(double r, double x) {
    // x ≡ x or y or z, the variable we are differentiating w.r.t.
    return m_1sNormalization * (-m_alpha) * x * exp(-m_alpha*r) / r;
}

double HydrogenOrbital::computeDerivative2s(double r, double x) {
    return m_2sNormalization * 0.25*m_alpha * (m_alpha * r - 4) * x * exp(-0.5* m_alpha*r) / r;
}

double HydrogenOrbital::computeDerivative2px(double  r,
                                             double  x,
                                             double  y,
                                             double  z,
                                             int     i) {
    const double N = m_2pNormalization;
    // i determines which variable we are differentiating w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return -N*(m_alpha*x*x - 2*r) * exp(-0.5*m_alpha*r) / (2*r);
    } else if (i==1) { // y derivative
        return -N*m_alpha*x*y * exp(-0.5*m_alpha*r) / (2*r);
    } else { // z derivative
        return -N*m_alpha*x*z * exp(-0.5*m_alpha*r) / (2*r);
    }
}

double HydrogenOrbital::computeDerivative2py(double r,
                                             double x,
                                             double y,
                                             double z,
                                             int    i) {
    // i determines which variable we are differentiating w.r.t., x, y, or z.
    const double N = m_2pNormalization;
    if (i==0) { // x derivative
        return -N*m_alpha*x*y * exp(-0.5*m_alpha*r) / (2*r);
    } else if (i==1) { // y derivative
        return -N*(m_alpha*y*y - 2*r) * exp(-0.5*m_alpha*r) / (2*r);
    } else { // z derivative
        return -N*m_alpha*y*z * exp(-0.5*m_alpha*r) / (2*r);
    }
}

double HydrogenOrbital::computeDerivative2pz(double r,
                                             double x,
                                             double y,
                                             double z,
                                             int i) {
    const double N = m_2pNormalization;
    // i determines which variable we are differentiating w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return -N*m_alpha*x*z * exp(-0.5*m_alpha*r) / (2*r);
    } else if (i==1) { // y derivative
        return -N*m_alpha*y*z * exp(-0.5*m_alpha*r) / (2*r);
    } else { // z derivative
        return -N*(m_alpha*z*z - 2*r) * exp(-0.5*m_alpha*r) / (2*r);
    }
}

double HydrogenOrbital::computeLaplacian1s(double r) {
    return m_1sNormalization * (m_alpha2 - 2 * m_alpha / r) * exp(-m_alpha * r);
}

double HydrogenOrbital::computeLaplacian2s(double r) {
    return m_2sNormalization * (1.25 * m_alpha2 - 2 * m_alpha / r - 0.125 * m_alpha2*m_alpha * r) * exp(-0.5* m_alpha * r);
}

double HydrogenOrbital::computeLaplacian2p(double r, double x) {
    // x ≡ x or y or z, for 2px, 2py, and 2pz, respectively.
    return m_2pNormalization * m_alpha*x * (m_alpha*r - 8.) * exp(-0.5*m_alpha*r) / (4*r);
}



























