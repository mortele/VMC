#include "heliumwithjastrow.h"
#include "WaveFunctions/wavefunction.h"
#include "system.h"
#include "electron.h"

using arma::vec;
using arma::mat;
using arma::norm;
using arma::zeros;
using arma::dot;
using std::exp;

HeliumWithJastrow::HeliumWithJastrow(System* system,
                                     double  alpha,
                                     double  beta,
                                     bool    useNumericalDerivatives) :
        HeliumWaveFunction(system, alpha, useNumericalDerivatives) {
    m_beta = beta;
}

double HeliumWithJastrow::evaluateWaveFunction() {
    const double base    = HeliumWaveFunction::evaluateWaveFunction();
    const double r12 = norm(m_system->getElectrons().at(0)->getPosition() -
                            m_system->getElectrons().at(1)->getPosition());
    const double jastrow = exp(0.5 * r12 / (1 + m_beta * r12));
    return base * jastrow;
}

double HeliumWithJastrow::evaluateLaplacian() {
    if (! m_useNumericalDerivatives) {
        // TODO: FIX THIS.
        const vec position1 = m_system->getElectrons().at(0)->getPosition();
        const vec position2 = m_system->getElectrons().at(1)->getPosition();
        const double r1  = norm(position1);
        const double r2  = norm(position2);
        const double r12        = norm(position1 - position2);
        const double r12Inverse = 1.0 / r12;
        const double r1Dotr2    = dot(position1, position2);
        const double r1Inverse  = 1.0 / r1;
        const double r2Inverse  = 1.0 / r2;

        const double a = m_alpha;
        const double b = m_beta;
        const double br12Plus1 = 1 + b*r12;

        const double A = 2*a*a - 2*a*(1/r1 + 1/r2);
        const double B = 0.5 / (br12Plus1*br12Plus1);
        const double C = -2*a/r12 * (r1+r2) * (1 - r1Dotr2/(r1*r2));
        const double D = 4/r12 - 4*b/br12Plus1 + 2/(br12Plus1*br12Plus1);
        const double L = A+B*(C+D);
        return L;
    } else {
        return WaveFunction::evaluateLaplacian();
    }
}








