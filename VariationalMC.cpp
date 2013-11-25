#include <cstdlib>
#include <cmath>

#include <iostream>

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
    int i;
    int j;
    static double local_energy;

    this->psi = trial_psi;
    this->var_params = var_params;
    this->N = n_points;

    energy = 0;
    variance = 0;
    accepts = 0;

    (*psi).set_params(var_params);

    set_seed(getUnixTime()*1000);

    for( i = 0; i < n_points; i++ )
    {
        for( j = 0; j < dimension; j++ )
        {
            cout << R[j] << " ";
        }
        cout << endl;
        get_newpos(R,R_new);

        accept_test(R,R_new);

        local_energy = get_localenergy(R);

        energy += local_energy;
        variance += local_energy*local_energy;
    }

    mean_energy = energy/n_points;
    std_energy = sqrt((variance/n_points -
                       mean_energy*mean_energy)/n_points);

    return mean_energy;
}

double VariationalMC::get_std()
{
    return std_energy;
}

double VariationalMC::get_acceptance_rate()
{
    return (double) accepts/N;
}


double VariationalMC::get_localenergy(double *R)
{
    return (*H)(psi, R)/(*psi)(R);
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

    if( rel_prob > 1 )
    {
        accept = true;
    }
    else
    {
        random_num = (double)rand() /  RAND_MAX;

        if( random_num <= rel_prob )
        {
            accept = true;
        }
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

