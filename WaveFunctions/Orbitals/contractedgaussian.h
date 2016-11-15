#pragma once
#include <vector>

class ContractedGaussian {
private:
    int                                     m_numberOfPrimitives;
    double                                  m_x;
    double                                  m_y;
    double                                  m_z;
    std::vector<class PrimitiveGaussian*>   m_primitives;

public:
    ContractedGaussian(int numberOfPrimitives, double x, double y, double z);
    void addPrimitive(class PrimitiveGaussian* primitive);
    double operator()(double x, double y, double z);
};

