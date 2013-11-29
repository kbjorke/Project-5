#include <cstdlib>
#include <cmath>

#include <iostream>

#include "mpi.h"
#include "VariationalMC.h"
#include "UnixTime.h"

using namespace std;


VariationalMC::VariationalMC(int dimension)
{
    this->dimension = dimension;

    R = new double[dimension];
    R_new = new double[dimension];
}

void VariationalMC::initialize(Hamiltonian *H,
                               double delta, double *R_init)
{
    int i;

    this->H = H;
    this->delta = delta;

    for( i = 0; i < dimension; i++ )
    {
        R[i] = R_init[i];
    }
}

double VariationalMC::operator()(Function *trial_psi,
                                 double *var_params, int n_points)
{
    int i, j; //, numprocs, my_rank;
    static double loc_en, local_energy, local_variance;

    this->psi = trial_psi;
    this->var_params = var_params;
    this->N = n_points;

    local_energy = 0;
    local_variance = 0;

    energy = 0;
    variance = 0;
    accepts = 0;

    (*psi).set_params(var_params);

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    set_seed(getUnixTime()*1000+my_rank);

    j = 0;
    for( i = my_rank; i < n_points; i += numprocs )
    {
        get_newpos(R,R_new);

        accept_test(R,R_new);

        /*
        for( k = 0; k < dimension; k++ )
        {
            cout << R[k] << " ";
        }
        cout << endl;
        cout << get_localenergy(R) << endl;
        */

        loc_en = get_localenergy(R);
        local_energy += loc_en;
        local_variance += loc_en*loc_en;
    }

    MPI_Reduce(&local_energy, &energy, 1, MPI_DOUBLE, MPI_SUM,
           0, MPI_COMM_WORLD);
    MPI_Reduce(&local_variance, &variance, 1, MPI_DOUBLE, MPI_SUM,
           0, MPI_COMM_WORLD);

    if( my_rank == 0 )
    {
        mean_energy = energy/n_points;
        std_energy = sqrt((variance/n_points -
                           mean_energy*mean_energy)/n_points);

        return mean_energy;
    }
    MPI_Finalize();

}

double VariationalMC::get_std()
{
    return std_energy;
}

double VariationalMC::get_acceptance_rate()
{
    return (double) accepts/(N);
}


double VariationalMC::get_localenergy(double *R)
{
    return (*H)(psi, R) / (*psi)(R);
}

void VariationalMC::get_newpos(double *R, double *R_new)
{
    static int i;
    static double random_num;

    for( i = 0; i < dimension; i++ )
    {
        random_num = 2*((double)rand() /  RAND_MAX) - 1;
        R_new[i] = R[i] + delta*random_num;
    }
}

void VariationalMC::accept_test(double *R, double *R_new)
{
    static int i;
    static double rel_prob, random_num;
    static bool accept;

    accept = false;

    rel_prob = relative_probability(R, R_new);

    random_num = (double)rand() /  RAND_MAX;

    if( random_num <= rel_prob )
    {
        accept = true;
    }

    if( accept )
    {
        for( i = 0; i < dimension; i++ )
        {
            R[i] = R_new[i];
        }

        accepts += 1;
    }
}

double VariationalMC::relative_probability(double *R, double *R_new)
{
    return ( (*psi)(R_new)*(*psi)(R_new) ) /
            ( (*psi)(R)*(*psi)(R) );
}

void VariationalMC::set_seed(double seed)
{
    srand(seed);
}
