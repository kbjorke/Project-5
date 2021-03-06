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

void MetropolisQM::initialize(Function *psi, double delta)
{
    int i;

    this->psi = psi;
    this->delta = delta;
}

double MetropolisQM::operator()(Observable *O, int n_points)
{
    static long int time_int;
    static int i, therm;
    static double loc_en, local_energy, local_variance;

    this->O = O;
    this->N = n_points;

    local_energy = 0;
    local_variance = 0;
    accepts = 0;

    energy = 0;
    variance = 0;
    total_accepts = 0;

    time_int = (long int) (getUnixTime()*10000 + my_rank);

    set_seed(time_int);

    // Set random starting points
    for( i = 0; i < 6; i++ )
    {
        R[i] = 2*ran2(&idum) - 1;
    }

    for( i = my_rank; i < n_points+thermalization; i += numprocs )
    {
        get_newpos(R,R_new);

        accept_test(R,R_new);

        if( i >= thermalization )
        {
            loc_en = get_localeigenvalue(R);
            local_energy += loc_en;
            local_variance += loc_en*loc_en;
        }
    }

    MPI_Reduce(&local_energy, &energy, 1, MPI_DOUBLE, MPI_SUM,
           0, MPI_COMM_WORLD);
    MPI_Reduce(&local_variance, &variance, 1, MPI_DOUBLE, MPI_SUM,
           0, MPI_COMM_WORLD);
    MPI_Reduce(&accepts, &total_accepts, 1, MPI_INT, MPI_SUM,
           0, MPI_COMM_WORLD);

    if( my_rank == 0 )
    {
        mean_energy = energy/n_points;
        std_energy = sqrt((variance/n_points -
                           mean_energy*mean_energy)/n_points);

        return mean_energy;
    }
    else
    {
        return 0;
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


double MetropolisQM::get_localeigenvalue(double *R)
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
    static long int seed_;

    seed_ = seed - ((long int) seed/(int) 1e7)*(int) 1e7;
    idum = -seed_;
}

void MetropolisQM::set_thermalization(int thermalization)
{
    this->thermalization = thermalization;
}
