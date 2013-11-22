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

#endif // PROBLEM_DEFINITIONS_H
