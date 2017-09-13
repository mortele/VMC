#include "slatertypeorbital.h"
#include <cmath>
#include "WaveFunctions/Orbitals/orbital.h"

SlaterTypeOrbital::SlaterTypeOrbital(double alpha) :
        Orbital() {

    double pi = acos(-1);

    m_alpha  = alpha;
    m_alpha2 = m_alpha*m_alpha;
    double a3 = pow(m_alpha,3);
    double a5 = pow(m_alpha,5);
    double a7 = pow(m_alpha,7);
    m_1sNormalization = sqrt(a3/pi);
    m_2sNormalization = sqrt(a5/(3*pi));
    m_2pNormalization = sqrt(2*a7/(15*pi));
}

double SlaterTypeOrbital::evaluate1s(double r) {
    return m_1sNormalization * exp(-m_alpha * r);
}

double SlaterTypeOrbital::evaluate2s(double r) {
    return m_2sNormalization * r * exp(-m_alpha * r);
}

double SlaterTypeOrbital::evaluate2p(double r, double x) {
    return m_2pNormalization * r * x * exp(-m_alpha * r);
}

double SlaterTypeOrbital::computeDerivative1s(double r, double x) {
    return evaluate1s(r) * x/r;
}

double SlaterTypeOrbital::computeDerivative2s(double r, double x) {
    return m_2sNormalization * x * exp(-m_alpha * r) * (1/r - m_alpha);
}

double SlaterTypeOrbital::computeDerivative2px(double r,
                                               double x,
                                               double y,
                                               double z,
                                               int i) {
    const double N = m_2pNormalization;
    // i determines which variable we are differentiating w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return N*exp(-m_alpha * r) / r * (r*r + x*x*(1 - r*m_alpha));
    } else if (i==1) { // y derivative
        return N*exp(-m_alpha * r) * x * y / r * (1 - r*m_alpha);
    } else { // z derivative
        return N*exp(-m_alpha * r) * x * z / r * (1 - r*m_alpha);
    }
}

double SlaterTypeOrbital::computeDerivative2py(double r,
                                               double x,
                                               double y,
                                               double z,
                                               int i) {
    const double N = m_2pNormalization;
    // i determines which variable we are differentiating w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return N*exp(-m_alpha * r) * x * y / r * (1 - r*m_alpha);
    } else if (i==1) { // y derivative
        return N*exp(-m_alpha * r) / r * (r*r + y*y*(1 - r*m_alpha));
    } else { // z derivative
        return N*exp(-m_alpha * r) * y * z / r * (1 - r*m_alpha);
    }
}

double SlaterTypeOrbital::computeDerivative2pz(double r,
                                               double x,
                                               double y,
                                               double z,
                                               int i) {
    const double N = m_2pNormalization;
    // i determines which variable we are differentiating w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return N*exp(-m_alpha * r) * x * z / r * (1 - r*m_alpha);
    } else if (i==1) { // y derivative
        return N*exp(-m_alpha * r) * y * z / r * (1 - r*m_alpha);
    } else { // z derivative
        return N*exp(-m_alpha * r) / r * (r*r + z*z*(1 - r*m_alpha));
    }
}

double SlaterTypeOrbital::computeLaplacian1s(double r) {
    // N = √a^3 / √pi
    return evaluate1s(r) * m_alpha * (r*m_alpha - 2) / r;
}

double SlaterTypeOrbital::computeLaplacian2s(double r) {
    // N = √a^5 / √(3 pi)
    return evaluate2s(r) / (r*r) * (2 - 4*r*m_alpha + r*r*m_alpha2);
}

double SlaterTypeOrbital::computeLaplacian2p(double r, double x) {
    // N = √(2 a^7) / √(15 pi)
    return m_2pNormalization * exp(-m_alpha * r) * x / r  * (4 - 6*r*m_alpha + r*r*m_alpha2);
}

double SlaterTypeOrbital::evaluate(double x,
                                   double y,
                                   double z,
                                   int index,
                                   int) {
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
        cout << "0 Invalid slater type orbital index: " << index << endl;
        exit(1);
    }
}

double SlaterTypeOrbital::computeDerivativeX(double x,
                                             double y,
                                             double z,
                                             int index) {
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
        cout << "1 Invalid slater type orbital index: " << index << endl;
        exit(1);
    }
}

double SlaterTypeOrbital::computeDerivativeY(double x,
                                             double y,
                                             double z,
                                             int index) {
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
        cout << "2 Invalid slater type orbital index: " << index << endl;
        exit(1);
    }
}

double SlaterTypeOrbital::computeDerivativeZ(double x,
                                             double y,
                                             double z,
                                             int index) {
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
        cout << "3 Invalid slater type orbital index: " << index << endl;
        exit(1);
    }
}

double SlaterTypeOrbital::computeDerivative(double x,
                                            double y,
                                            double z,
                                            int index,
                                            int dimension,
                                            int ) {
    if (dimension==0) {
        return computeDerivativeX(x,y,z,index);
    } else if (dimension==1) {
        return computeDerivativeY(x,y,z,index);
    } else {
        return computeDerivativeZ(x,y,z,index);
    }
}

double SlaterTypeOrbital::computeLaplacian(double x,
                                           double y,
                                           double z,
                                           int index,
                                           int ) {
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
        cout << "4 Invalid slater type orbital index: " << index << endl;
        exit(1);
    }
}
