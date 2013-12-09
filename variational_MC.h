#ifndef VARIATIONAL_MC_H
#define VARIATIONAL_MC_H

#include <sstream>

#include "Function.h"
#include "Observable.h"

using namespace std;

/* Function to preform Variational MonteCarlo for a given system
 * and a given trial wave function. The variation will be preform
 * on an intervall given by the user.
 * */
double variational_MC(Function *psi_trial, Observable *hamiltonian,
                      double *params, int num_params, double *var_params,
                      int ind_var, int num_var, double *best_var, int N,
                      string *id_params, bool output = false,
                      string trial_lable = "none");

#endif // VARIATIONAL_MC_H
