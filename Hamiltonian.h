#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "Function.h"
#include "Observable.h"

/* Class for implementation of Hamiltonian.
 *
 * Used to create a callable object to represent the Hamiltonian
 * to be used for Quantum Mechanical problems. Subclass of
 * Observable superclass.
 *
 * Methods:
 *          - contructor     : Used to create Hamiltonian object.
 *          - operator()     : Used to apply the Hamiltonian to a state or
 *                             wavefunction at a specific point.
 *          - set_potential  : Used to set a external potential to the
 *                             Hamiltonian.
 * */
class Hamiltonian : public Observable
{
protected:
    Function *V;

public:
    Hamiltonian();

    double operator()(Function *psi, double *args);
    void set_potential(Function *V);
};

#endif // HAMILTONIAN_H
