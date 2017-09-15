#pragma once
#include <cmath>
#include "WaveFunctions/Orbitals/orbital.h"

class HydrogenOrbital : public Orbital {
private:
    double  m_alpha  = 1.0;
    double  m_alpha2 = 1.0;
    double  m_1sNormalization = 1;
    double  m_2sNormalization = 1;
    double  m_2pNormalization = 1;

    double evaluate1s (double r);
    double evaluate2s (double r);
    double evaluate2p (double r, double x);

    double computeDerivative1s(double r, double x);
    double computeDerivative2s(double r, double x);
    double computeDerivative2px(double r, double x, double y, double z, int i);
    double computeDerivative2py(double r, double x, double y, double z, int i);
    double computeDerivative2pz(double r, double x, double y, double z, int i);

    double computeLaplacian1s(double r);
    double computeLaplacian2s(double r);
    double computeLaplacian2p(double r, double x);

public:
    HydrogenOrbital(double alpha);

    double evaluate(double x, double y, double z, int index, int spin=0);
    double computeDerivativeX(double x, double y, double z, int index);
    double computeDerivativeY(double x, double y, double z, int index);
    double computeDerivativeZ(double x, double y, double z, int index);
    double computeDerivative(double x, double y, double z, int index, int dimension, int spin=0);
    double computeLaplacian(double x, double y, double z, int index, int spin=0);
};
