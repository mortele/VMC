#pragma once
#include <vector>
#include <iostream>

class ContractedGaussian {
    friend class PrimitiveGaussian;

private:
    int                                     m_numberOfPrimitives;
    double                                  m_x;
    double                                  m_y;
    double                                  m_z;
    double                                  m_currentValue;
    std::vector<class PrimitiveGaussian*>   m_primitives;

public:
    ContractedGaussian(int numberOfPrimitives, double x, double y, double z);
    void addPrimitive(class PrimitiveGaussian* primitive);
    double operator()(double x, double y, double z);
    double evaluate(double x, double y, double z);

    double xDerivative(double x, double y, double z);
    double yDerivative(double x, double y, double z);
    double zDerivative(double x, double y, double z);

    double calculateLaplacian(double x, double y, double z);
    friend std::ostream& operator<<(std::ostream& stream, const ContractedGaussian& contracted);
};

