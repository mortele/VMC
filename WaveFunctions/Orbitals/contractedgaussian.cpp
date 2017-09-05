#include "contractedgaussian.h"
#include "WaveFunctions/Orbitals/primitivegaussian.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;
using std::setprecision;


ContractedGaussian::ContractedGaussian(int      numberOfPrimitives,
                                       double   x,
                                       double   y,
                                       double   z) {
    m_numberOfPrimitives = numberOfPrimitives;
    m_primitives.reserve(numberOfPrimitives);
    m_x = x;
    m_y = y;
    m_z = z;
}

void ContractedGaussian::addPrimitive(PrimitiveGaussian* primitive) {
    m_primitives.push_back(primitive);
}

double ContractedGaussian::operator()(double x, double y, double z) {
    double result = 0;
    for (PrimitiveGaussian* primitive : m_primitives) {
        result += (*primitive)(x - m_x, y - m_y, z - m_z);
    }
    m_currentValue = result;
    return result;
}

double ContractedGaussian::evaluate(double x, double y, double z) {
    return (*this)(x,y,z);
}

double ContractedGaussian::xDerivative(double x, double y, double z) {
    /*cout << x << " " << y << " " << z << endl;
    cout << "=======================" << endl;
    cout << std::sqrt(x*x+y*y+z*z) << endl;
    cout << "=======================" << endl;*/
    evaluate(x,y,z);
    double result = 0;
    for (int i = 0; i < m_numberOfPrimitives; i++) {
        PrimitiveGaussian*& primitive = m_primitives.at(i);
        result += primitive->xDerivative(x,y,z);// * primitive->getCurrentValue();

        /*printf("%10.6g (%d,%d,%d) exp(-%10.6g r²)   -->   %10g   --> dx: %10g\n",
               primitive->m_constant,
               primitive->m_i,
               primitive->m_j,
               primitive->m_k,
               primitive->m_alpha,
               primitive->getCurrentValue(),
               primitive->xDerivative(x,y,z));*/
    }

    double tmp = result;// /m_currentValue;
    //cout << "value: " << m_currentValue << endl;
    //cout << "total: " << tmp << endl;



    //exit(1);
    return tmp;
}

double ContractedGaussian::yDerivative(double x, double y, double z) {
    double result = 0;
    for (int i = 0; i < m_numberOfPrimitives; i++) {
        PrimitiveGaussian*& primitive = m_primitives.at(i);
        result += primitive->yDerivative(x,y,z);// * primitive->getCurrentValue();
    }
    return result;// /m_currentValue;
}

double ContractedGaussian::zDerivative(double x, double y, double z) {
    double result = 0;
    for (int i = 0; i < m_numberOfPrimitives; i++) {
        PrimitiveGaussian*& primitive = m_primitives.at(i);
        result += primitive->zDerivative(x,y,z);// * primitive->getCurrentValue();
    }
    return result;// /m_currentValue;
}

double ContractedGaussian::calculateLaplacian(double x, double y, double z) {
    /*cout << ".--.-.-.-.-.-.-.-.-.--." << endl;
    cout << ".--.-.-.-.-.-.-.-.-.--." << endl;
    cout << x << " " << y << " " << z << endl;
    cout << "=======================" << endl;
    cout << std::sqrt(x*x+y*y+z*z) << endl;
    cout << "=======================" << endl;*/

    evaluate(x,y,z);
    const double xA         = x - m_x;
    const double yA         = y - m_y;
    const double zA         = z - m_z;
    //double totalValue       = 0;
    double totalLaplacian   = 0;
    for (PrimitiveGaussian* primitive : m_primitives) {
        const double ddx     = primitive->xxDerivative(xA,yA,zA);
        const double ddy     = primitive->yyDerivative(xA,yA,zA);
        const double ddz     = primitive->zzDerivative(xA,yA,zA);
        totalLaplacian      += (ddx + ddy + ddz) * primitive->getCurrentValue();

        /*printf("%10.6g (%d,%d,%d) exp(-%10.6g r²)   -->   %10g   --> dx: %10g\n",
               primitive->m_constant,
               primitive->m_i,
               primitive->m_j,
               primitive->m_k,
               primitive->m_alpha,
               primitive->getCurrentValue(),
               primitive->xxDerivative(x,y,z));*/
    }
    //m_currentValue = totalValue;
    //cout << m_currentValue << endl;
    //cout << totalLaplacian << endl;//   / m_currentValue << endl;
    //exit(1);
    return totalLaplacian;// / m_currentValue;
}








