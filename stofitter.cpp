#include "stofitter.h"
#include "Optimization/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Optimization/interpolation.h"
#include "RandomNumberGenerator/random.h"

using namespace alglib;
using std::cout;
using std::endl;

int STOFitter::m_n;


void STOFitter::G(const real_1d_array &c,
                  const real_1d_array &x,
                  double              &function,
                  void                *ptr) {

    double r2 = x[0]*x[0];
    function = c[0] *exp(-c[1] * r2);
    for (int ii = 2; ii < 2*m_n; ii+=2) {
        function += c[ii] *exp(-c[ii+1] * r2);
    }
    //func *= r2;
}

void STOFitter::gradient(const real_1d_array &c,
                         const real_1d_array &x,
                         double              &function,
                         real_1d_array       &gradient,
                         void*               ptr) {

    double r2 = x[0]*x[0];
    function = c[0] *exp(-c[1] * r2);
    for (int ii = 2; ii < 2*m_n; ii+=2) {
        function += c[ii] *exp(-c[ii+1] * r2);
    }
    //func *= r2;

    gradient[0] = (         exp(-c[1] * r2));
    gradient[1] = (-r2*c[0]*exp(-c[1] * r2));
    for (int ii = 2; ii < 2*m_n; ii+=2) {
        gradient[ii]   = (          exp(-c[ii+1] * r2));
        gradient[ii+1] = (-r2*c[ii]*exp(-c[ii+1] * r2));
    }
}

void STOFitter::hessian(const real_1d_array &c,
                        const real_1d_array &x,
                        double              &function,
                        real_1d_array       &gradient,
                        real_2d_array       &hessian,
                        void*               ptr) {

    double r2 = x[0]*x[0];
    function = c[0] *exp(-c[1] * r2);
    for (int ii = 2; ii < 2*m_n; ii+=2) {
        function += c[ii] *exp(-c[ii+1] * r2);
    }
    //func *= r2;

    gradient[0] = (         exp(-c[1] * r2));
    gradient[1] = (-r2*c[0]*exp(-c[1] * r2));
    for (int ii = 2; ii < 2*m_n; ii+=2) {
        gradient[ii]   = (          exp(-c[ii+1] * r2));
        gradient[ii+1] = (-r2*c[ii]*exp(-c[ii+1] * r2));
    }

    hessian[0][0] = c[0]*r2*r2*exp(-c[1]*r2);
    hessian[0][1] = -r2*exp(-c[1]*r2);
    hessian[1][0] = -r2*exp(-c[1]*r2);
    hessian[1][1] = 0;
    for (int ii = 2; ii < 2*m_n; ii+=2) {
        hessian[ii][ii]     = 0;
        hessian[ii][ii+1]   = -r2*exp(-c[ii+1]*r2);
        hessian[ii+1][ii]   = -r2*exp(-c[ii+1]*r2);
        hessian[ii+1][ii+1] = c[ii]*r2*r2*exp(-c[ii+1]*r2);
    }

}

double STOFitter::s1(double r) {
    //return sqrt(pow(m_alpha,3)/m_pi) *
    return exp(-m_alpha * r*r);// * r*r;
}

STOFitter::STOFitter(double alpha) {
    m_alpha = alpha;
    m_pi = acos(-1);
}

void STOFitter::computeFit(int n) {

    n = 1;

    STOFitter::m_n = n;
    real_1d_array c;
    real_1d_array l;
    real_1d_array u;
    double* c_ = new double[2*n];
    double* u_ = new double[2*n];
    double* l_ = new double[2*n];
    for (int i = 0; i < 2*n; i+=2) {
        c_[i]   = Random::nextDouble(0.0, 1.0);
        c_[i+1] = Random::nextDouble(0.0, 1.0);
        l_[i]   = -1e4;
        l_[i+1] = -1e4;
        u_[i]   = 1e4;
        u_[i+1] = 1e4;
    }
    c.setcontent(2*n,c_);
    l.setcontent(2*n,l_);
    u.setcontent(2*n,u_);

    int N = 1000;
    double x0 = 0;
    double x1 = 2;
    double dx = (x1-x0)/(N-1);
    double* y_ = new double[N];
    double* x_ = new double[N];
    for (int i = 0; i < N; i++) {
        x_[i] = x0+i*dx;
        y_[i] = s1(x_[i]);

    }


    real_1d_array y;
    real_2d_array x;
    y.setcontent(N,y_);
    x.setcontent(N,1,y_);



    double epsx = 1e-10;
    ae_int_t maxits = 1e4;
    ae_int_t info;
    lsfitstate state;
    lsfitreport rep;
    double dxx = 1e-4;

    lsfitcreatefg(x, y, c, dxx, state);
    lsfitsetcond(state, epsx, maxits);
    lsfitsetbc(state,l,u);
    alglib::lsfitfit(state, G, gradient);
    lsfitresults(state, info, c, rep);
    printf("%g\n", rep.r2);
    printf("%d\n", int(info));
    printf("%s\n", c.tostring(6).c_str());
}
