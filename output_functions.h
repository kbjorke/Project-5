#ifndef OUTPUT_FUNCTIONS_H
#define OUTPUT_FUNCTIONS_H

#include <sstream>

#include "Function.h"
#include "Observable.h"

using namespace std;

void test_integrators();
void output_method(string method);
/*
double variational_MC(Function *psi_trial, Observable *hamiltonian,
                      double * R_init, double *params,
                      int num_params, double *var_params, int ind_var,
                      int num_var, double *best_var, string *id_params,
                      int N);
                      */
void output_VMC_data_header(fstream *outfile, int N, int ind_var,
                            double *params, int num_params,
                            string *id_params);
void output_VMC_data(fstream *outfile, double *params, int ind_var,
                     double energy, double stdev, double acceptance_rate,
                     double delta_t);
void output_VMC_data_end(fstream *outfile, double best_energy,
                         string *id_params, int ind_best,
                         double best_var);

#endif // OUTPUT_FUNCTIONS_H
