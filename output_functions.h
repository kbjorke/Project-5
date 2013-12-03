#include <sstream>

#include "Function.h"
#include "Observable.h"

using namespace std;

void test_integrators();
void output_method(string method);
double variational_MC(Function *psi_trial, Observable *hamiltonian,
                      double * R_init, double *params,
                      int num_params, double *var_params, int ind_var,
                      int num_var, double *best_var, string *id_params,
                      int N);
void output_VMC_data_header();
void output_VMC_data();



