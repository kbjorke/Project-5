#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include "Function.h"

class Observable
{
public:
    Observable(){}
    virtual double operator()(Function *psi, double *args){}
    virtual void set_potential(Function *V){}
    virtual void set_params(double* params){}

    ~Observable(){}
};

#endif // OBSERVABLE_H
