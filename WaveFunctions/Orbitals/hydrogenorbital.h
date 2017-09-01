#pragma once
#include <cmath>

class HydrogenOrbital {
private:
    double  m_alpha  = 1.0;
    double  m_alpha2 = 1.0;
    int     m_index  = -1;

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

    double operator()(double x, double y, double z, int index);
    double evaluate(double x, double y, double z, int index);
    double computeDerivativeX(double x, double y, double z, int index);
    double computeDerivativeY(double x, double y, double z, int index);
    double computeDerivativeZ(double x, double y, double z, int index);
    double computeLaplacian(double x, double y, double z, int index);
};
