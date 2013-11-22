#include "Hamiltonian.h"

Hamiltonian::Hamiltonian()
{
}

double Hamiltonian::operator()(Function *psi, double *args)
{
    static double K;

    K = - 0.5*psi->derivative2(args);

    return K + (*V)(args);
}

void Hamiltonian::set_potential(Function *V)
{
    this->V = V;
}
