#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "Function.h"

class Hamiltonian
{
protected:
    Function *V;

public:
    Hamiltonian();

    double operator()(Function *psi, double *args);
    void set_potential(Function *V);

    virtual ~Hamiltonian(){}
};

#endif // HAMILTONIAN_H
