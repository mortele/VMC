#include "slatertypeorbital.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "WaveFunctions/Orbitals/orbital.h"

using std::cout;
using std::endl;
using std::setprecision;

SlaterTypeOrbital::SlaterTypeOrbital(double alpha) :
        Orbital() {

    double pi = acos(-1.);
    m_alpha  = alpha;
    m_alpha2 = m_alpha*m_alpha;
    double a3 = pow(m_alpha,3);
    double a5 = pow(m_alpha,5);
    double a7 = pow(m_alpha,7);
    m_1sNormalization = sqrt(a3/pi);
    m_2sNormalization = (1/4.)*sqrt(a5/(6*pi));
    m_2pNormalization = (1/8.)*sqrt(a7/(15*pi));
}

inline double SlaterTypeOrbital::evaluate1s(double r) {
    return m_1sNormalization * exp(-m_alpha * r);
}

inline double SlaterTypeOrbital::evaluate2s(double r) {
    return m_2sNormalization * r * exp(-0.5 * m_alpha * r);
}

inline double SlaterTypeOrbital::evaluate2p(double r, double x) {
    const double a = 0.5*m_alpha;
    return m_2pNormalization * r * x * exp(-a * r);
}

inline double SlaterTypeOrbital::computeDerivative1s(double r, double x) {
    return -evaluate1s(r) * m_alpha * x / r ;
}

inline double SlaterTypeOrbital::computeDerivative2s(double r, double x) {
    //return m_2sNormalization * x * exp(-0.5 * m_alpha * r) * (m_alpha -2/r) / 2.;
    return - sqrt(pow(m_alpha,5)/(6*acos(-1.))) * exp(-0.5*m_alpha*r) * x * (-2 + r*m_alpha) / (8 * r);
}

inline double SlaterTypeOrbital::computeDerivative2px(double r,
                                               double x,
                                               double y,
                                               double z,
                                               int i) {
    const double a = 0.5*m_alpha;
    const double N = m_2pNormalization;
    // i determines which variable we are differentiating w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return N*exp(-a * r) / r * (r*r + x*x*(1 - r*a));
    } else if (i==1) { // y derivative
        return N*exp(-a * r) * x * y / r * (1 - r*a);
    } else { // z derivative
        return N*exp(-a * r) * x * z / r * (1 - r*a);
    }
}

inline double SlaterTypeOrbital::computeDerivative2py(double r,
                                               double x,
                                               double y,
                                               double z,
                                               int i) {
    const double a = 0.5*m_alpha;
    const double N = m_2pNormalization;
    // i determines which variable we are differentiating w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return N*exp(-a * r) * x * y / r * (1 - r*a);
    } else if (i==1) { // y derivative
        return N*exp(-a * r) / r * (r*r + y*y*(1 - r*a));
    } else { // z derivative
        return N*exp(-a * r) * y * z / r * (1 - r*a);
    }
}

inline double SlaterTypeOrbital::computeDerivative2pz(double r,
                                               double x,
                                               double y,
                                               double z,
                                               int i) {
    const double a = 0.5*m_alpha;
    const double N = m_2pNormalization;
    // i determines which variable we are differentiating w.r.t., x, y, or z.
    if (i==0) { // x derivative
        return N*exp(-a * r) * x * z / r * (1 - r*a);
    } else if (i==1) { // y derivative
        return N*exp(-a * r) * y * z / r * (1 - r*a);
    } else { // z derivative
        return N*exp(-a * r) / r * (r*r + z*z*(1 - r*a));
    }
}

inline double SlaterTypeOrbital::computeLaplacian1s(double r) {
    // N = √a^3 / √pi
    return m_1sNormalization * (m_alpha2 - 2 * m_alpha / r) * exp(-m_alpha * r);
}

inline double SlaterTypeOrbital::computeLaplacian2s(double r) {
    // N = √a^5 / 4√(6 pi)
    //return evaluate2s(r) / (r*r) * (2 - 4*r*m_alpha + r*r*m_alpha2);
    return m_2sNormalization / (4 * r) * exp(-0.5*m_alpha*r) * (8 - 8*r*m_alpha + r*r*m_alpha2);
}

inline double SlaterTypeOrbital::computeLaplacian2p(double r, double x) {
    // N = √(2 a^7) / √(15 pi)
    const double a = 0.5*m_alpha;
    return m_2pNormalization * exp(-a * r) * x / r  * (4 - 6*r*a + r*r*a*a);
}

inline double SlaterTypeOrbital::evaluate(double x,
                                   double y,
                                   double z,
                                   int index,
                                   int spin) {
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

inline double SlaterTypeOrbital::computeDerivativeX(double x,
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

inline double SlaterTypeOrbital::computeDerivativeY(double x,
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

inline double SlaterTypeOrbital::computeDerivativeZ(double x,
                                             double y,
                                             double z,
                                             int index) {
    //const double r = sqrt(x*x + y*y + z*z);
    if (index==0) {
        return computeDerivative1s(sqrt(x*x + y*y + z*z), z);
    } else if (index==1) {
        return computeDerivative2s(sqrt(x*x + y*y + z*z), z);
    } else if (index==2) {
        return computeDerivative2px(sqrt(x*x + y*y + z*z), x, y, z, 2);
    } else if (index==3) {
        return computeDerivative2py(sqrt(x*x + y*y + z*z), x, y, z, 2);
    } else if (index==4) {
        return computeDerivative2pz(sqrt(x*x + y*y + z*z), x, y, z, 2);
    } else {
        //cout << "3 Invalid slater type orbital index: " << index << endl;
        return 0;
    }
}

inline double SlaterTypeOrbital::computeDerivative(double x,
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

inline double SlaterTypeOrbital::computeLaplacian(double x,
                                           double y,
                                           double z,
                                           int index,
                                           int ) {
    //const double r = sqrt(x*x + y*y + z*z);
    if (index==0) {
        return computeLaplacian1s(sqrt(x*x + y*y + z*z));
    } else if (index==1) {
        return computeLaplacian2s(sqrt(x*x + y*y + z*z));
    } else if (index==2) {
        return computeLaplacian2p(sqrt(x*x + y*y + z*z), x);
    } else if (index==3) {
        return computeLaplacian2p(sqrt(x*x + y*y + z*z), y);
    } else if (index==4) {
        return computeLaplacian2p(sqrt(x*x + y*y + z*z), z);
    } else {
        return 0;
    }
}
