#include <cstdlib>

#include "VariationalMC.h"


VariationalMC::VariationalMC()
{
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
        random_num = (double)rand() /  RAND_MAX;
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

void VariationalMC::set_variational_params(double *var_params)
{
    this->var_params = var_params;
}
