#ifndef PROBLEM_DEFINITIONS_H
#define PROBLEM_DEFINITIONS_H

#include "Function.h"

class Integrand: public Function
{
    private:
    double alpha, mu, rho;
    double dist[3];

    public:
    Integrand(int dimension);
    double set_params(double *params);
    double operator()(double *r);
};

class Psi_T1: public Function
{
    private:
    double alpha, mu;

    public:
    Psi_T1(int dimension);
    double set_params(double *params);
    double operator()(double *r);
};

class Psi_T2: public Function
{
    private:
    double alpha, beta, mu, rho, exponent;
    double dist[3];

    public:
    Psi_T2(int dimension);
    double set_params(double *params);
    double operator()(double *r);
};

class V_TOT: public Function
{
    private:
    double mu, rho;
    double dist[3];

    public:
    V_TOT(int dimension);
    double operator()(double *r);
};

class V_HO: public Function
{
    private:
    double mu;

    public:
    V_HO(int dimension);
    double operator()(double *r);
};

#endif // PROBLEM_DEFINITIONS_H
