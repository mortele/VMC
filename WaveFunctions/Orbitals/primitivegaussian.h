#pragma once

class PrimitiveGaussian {
    friend class ContractedGaussian;

private:
    int     m_i;
    int     m_j;
    int     m_k;
    double  m_Ax;
    double  m_Ay;
    double  m_Az;
    double  m_alpha;
    double  m_constant;

    double m_currentValue;

    double pow(double a, int n);

public:
    PrimitiveGaussian(int    i,
                      int    j,
                      int    k,
                      double Ax,
                      double Ay,
                      double Az,
                      double alpha,
                      double constant);
    double operator()(double x, double y, double z);
    double getCurrentValue() { return m_currentValue; }

    double xDerivative(double x, double y, double z);
    double yDerivative(double x, double y, double z);
    double zDerivative(double x, double y, double z);

    double xxDerivative(double x, double y, double z);
    double yyDerivative(double x, double y, double z);
    double zzDerivative(double x, double y, double z);
};
