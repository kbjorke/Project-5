#ifndef PROBLEM_DEFINITIONS_H
#define PROBLEM_DEFINITIONS_H

#include "Function.h"
#include "Observable.h"


/* Contains classes which describe certain intgrands,
 * wavefunctions, Hamiltonians or wavefunctions for the project.
 * */
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

class Integrand2: public Function
{
    private:
    double rho;
    double dist[3];

    public:
    Integrand2(int dimension);
    double operator()(double *r);
};

class Psi_T1: public Function
{
    private:
    double alpha, alpha2, mu;

    public:
    Psi_T1(int dimension);
    double set_params(double *params);
    double operator()(double *r);
    double derivative2(double *r);
};

class Psi_T2: public Function
{
    private:
    double alpha, alpha2, beta, beta2, mu, rho, exponent;
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

class PotentialEnergy : public Observable
{
    private:
    double rho;
    double dist[3];

    public:
    PotentialEnergy(){}
    double operator()(Function *psi, double *r);
};

#endif // PROBLEM_DEFINITIONS_H
