#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include "Function.h"

/* Superclass for implementation of an Observable.
 *
 * Used as a superclass for a class that represents an observable
 * to be used for Quantum Mechanical problems.
 *
 * Methods:
 *          - contructor     : Used to create Observable object.
 *          - operator()     : Used to apply the Observable operator to a
 *                             state or wavefunction at a specific point,
 *                             given by *args.
 *          - set_potential  : Used to add a potential to the to the
 *                             observable.
 *          - set_params     : Used to set the parameters of the observable.
 * */
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
