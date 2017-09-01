#pragma once
#include <cmath>

class HydrogenOrbital {
private:
    double m_alpha  = 1.0;
    double m_alpha2 = 1.0;

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
    double computeDoubleDerivative2p(double r);

public:
    HydrogenOrbital(double alpha);
};
