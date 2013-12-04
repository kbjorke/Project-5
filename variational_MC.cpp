#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <ctime>

#include <iostream>

#include "mpi.h"
#include "globals.h"
#include "output_functions.h"
#include "MetropolisQM.h"
#include "UnixTime.h"

using namespace std;
double variational_MC(Function *psi_trial, Observable *hamiltonian,
                      double *R_init, double *params, int num_params,
                      double *var_params, int ind_var, int num_var,
                      double *best_var, int N, string *id_params,
                      bool output)
{
    int i, ind_best;
    double delta, stdev, t0, t1, delta_t, acceptance_rate;
    double energy, best_energy;

    fstream outfile;

    if( my_rank == 0 && output )
    {
        output_VMC_data_header(&outfile, N, ind_var, params,
                               num_params, id_params);
    }

    MetropolisQM metroQM(6);

    delta = 0.6;

    best_energy = 1e10;

    for( i = 0; i < num_var; i++ )
    {
        params[ind_var] = var_params[i];

        (*psi_trial).set_params(params);

        metroQM.initialize(psi_trial, delta, R_init);

        if( my_rank == 0 )
        {
            t0 = getUnixTime();
        }

        energy = metroQM(hamiltonian, N);

        if( my_rank == 0 )
        {
            t1 = getUnixTime();

            stdev = metroQM.get_std();
            acceptance_rate = metroQM.get_acceptance_rate();
            delta_t = t1 - t0;


            if( energy < best_energy )
            {
                best_energy = energy;
                ind_best = ind_var;
                *best_var = params[ind_var];
            }
        }

        if( my_rank == 0 && outfile )
        {
            output_VMC_data(&outfile, params, ind_var,
                            energy, stdev, acceptance_rate,
                            delta_t);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    if( my_rank == 0 & output )
    {
        output_VMC_data_end(&outfile, best_energy, id_params,
                            ind_best, *best_var);
    }

    return best_energy;
}
