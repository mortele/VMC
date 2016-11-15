#pragma once

class PrimitiveGaussian {
    friend class ContractedGaussian;

private:
    int     m_x;
    int     m_y;
    int     m_z;
    double  m_alpha;
    double  m_constant;

    double pow(double a, int n);

public:
    PrimitiveGaussian(int    x,
                      int    y,
                      int    z,
                      double alpha,
                      double constant);
    double operator()(double x, double y, double z);

    double xDerivative(double x, double y, double z);
    double yDerivative(double x, double y, double z);
    double zDerivative(double x, double y, double z);
};
