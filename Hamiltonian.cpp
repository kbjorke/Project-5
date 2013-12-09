#include "Hamiltonian.h"

using namespace std;

Hamiltonian::Hamiltonian(){}

/* Method operator().
 *
 * Used to apply the Hamiltonian to a state or wavefunction
 * at a given position/arguments.
 *
 * Input:
 *          *psi  : Pointer to a function object containing the
 *                  wavefunction to apply the Hamiltonian to.
 *
 *          *args : Pointer or array containing the parameters
 *                  where the Hamiltonian and wavefunction will
 *                  interact.
 * */
double Hamiltonian::operator()(Function *psi, double *args)
{
    static double K;

    K = - 0.5*(*psi).derivative2(args);

    return K + (*V)(args)*(*psi)(args);
}

/* Method set_potential.
 *
 * Used to set the potential for the Hamiltonian.
 *
 * Input:
 *          *V : Pointer to a function object containing the
 *               potential to be added.
 *
 * NOTE: Might add implementation or seperate method to add a
 *       new/extra potentials to the Hamiltonian.
 * */
void Hamiltonian::set_potential(Function *V)
{
    this->V = V;
}
