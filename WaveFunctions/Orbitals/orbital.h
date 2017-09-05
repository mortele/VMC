#pragma once

class Orbital {
public:
    Orbital() {}
    double operator()(double x, double y ,double z, int index) { return evaluate(x,y,z,index); }
    virtual double evaluate(double x,double y,double z,int index) = 0;
    virtual double computeDerivative(double x,double y,double z,int index,int dimension) = 0;
    virtual double computeLaplacian(double x,double y,double z,int index) = 0;
};

