#ifndef VARIATIONAL_MC_H
#define VARIATIONAL_MC_H

#include <sstream>

#include "Function.h"
#include "Observable.h"

using namespace std;

double variational_MC(Function *psi_trial, Observable *hamiltonian,
                      double * R_init, double *params,
                      int num_params, double *var_params, int ind_var,
                      int num_var, double *best_var, int N,
                      string *id_params, bool output = false);

#endif // VARIATIONAL_MC_H
