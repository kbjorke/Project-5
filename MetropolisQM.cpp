#include <cstdlib>
#include <cmath>

#include <iostream>

#include "globals.h"
#include "mpi.h"
#include "lib.h"
#include "MetropolisQM.h"
#include "UnixTime.h"

using namespace std;


MetropolisQM::MetropolisQM(int dimension)
{
    this->dimension = dimension;

    R = new double[dimension];
    R_new = new double[dimension];

    thermalization = 10; //default thermalization
}

void MetropolisQM::initialize(Function *psi,
                              double delta, double *R_init)
{
    int i;

    this->psi = psi;
    this->delta = delta;

    for( i = 0; i < dimension; i++ )
    {
        R[i] = R_init[i];
    }
}

double MetropolisQM::operator()(Observable *O, int n_points)
{
    static int i, therm;
    static double loc_en, local_energy, local_variance;

    this->O = O;
    this->N = n_points;

    local_energy = 0;
    local_variance = 0;

    energy = 0;
    variance = 0;
    accepts = 0;

    set_seed((long int) (getUnixTime()*100 + my_rank));

    therm = thermalization;

    for( i = my_rank; i < n_points+therm; i += numprocs )
    {
        get_newpos(R,R_new);

        accept_test(R,R_new);

        if( i >= therm )
        {
            loc_en = get_localenergy(R);
            local_energy += loc_en;
            local_variance += loc_en*loc_en;
        }
    }

    MPI_Reduce(&local_energy, &energy, 1, MPI_DOUBLE, MPI_SUM,
           0, MPI_COMM_WORLD);
    MPI_Reduce(&local_variance, &variance, 1, MPI_DOUBLE, MPI_SUM,
           0, MPI_COMM_WORLD);
    MPI_Reduce(&accepts, &total_accepts, 1, MPI_DOUBLE, MPI_SUM,
           0, MPI_COMM_WORLD);
    MPI_Reduce(&therm, &thermalization, 1, MPI_DOUBLE, MPI_SUM,
           0, MPI_COMM_WORLD);

    if( my_rank == 0 )
    {
        thermalization = thermalization/numprocs;
        mean_energy = energy/n_points;
        std_energy = sqrt((variance/n_points -
                           mean_energy*mean_energy)/n_points);

        return mean_energy;
    }
}

double MetropolisQM::get_std()
{
    return std_energy;
}

double MetropolisQM::get_acceptance_rate()
{
    return (double) total_accepts/(N + thermalization);
}


double MetropolisQM::get_localenergy(double *R)
{
    return (*O)(psi, R) / (*psi)(R);
}

void MetropolisQM::get_newpos(double *R, double *R_new)
{
    static int i;
    static double random_num;

    for( i = 0; i < dimension; i++ )
    {
        random_num = 2*ran2(&idum) - 1;
        //random_num = 2*((double)rand() /  RAND_MAX) - 1;
        R_new[i] = R[i] + delta*random_num;
    }
}

void MetropolisQM::accept_test(double *R, double *R_new)
{
    static int i;
    static double rel_prob, random_num;
    static bool accept;

    accept = false;

    rel_prob = relative_probability(R, R_new);

    random_num = ran2(&idum);
    //random_num = (double)rand() /  RAND_MAX;

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

double MetropolisQM::relative_probability(double *R, double *R_new)
{
    return ( (*psi)(R_new)*(*psi)(R_new) ) /
            ( (*psi)(R)*(*psi)(R) );
}

void MetropolisQM::set_seed(long int seed)
{
    idum = -seed;
    //srand(seed);
}

void MetropolisQM::set_thermalization(int thermalization)
{
    this->thermalization = thermalization;
}
