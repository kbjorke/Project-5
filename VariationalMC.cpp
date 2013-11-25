#include "VariationalMC.h"

VariationalMC::VariationalMC()
{
}

double VariationalMC::get_localenergy(double *R)
{
    return (*H)(psi, R)/(*psi)(R);
}


void get_newpos(double *R, double *R_new);
void accept_test(double *R, double *R_new);
double relative_probability(double *R, double *R_new);

void set_variational_params(double *var_params);
