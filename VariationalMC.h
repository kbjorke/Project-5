#ifndef VARIATIONALMC_H
#define VARIATIONALMC_H

#include "Function.h"
#include "Hamiltonian.h"

class VariationalMC
{
private:
    double *var_params;
    Hamiltonian *H;
    Function *psi;

    double get_localenergy(double *R);
    void get_newpos(double *R, double *R_new);
    void accept_test(double *R, double *R_new);
    double relative_probability(double *R, double *R_new);

public:
    VariationalMC();
    void set_variational_params(double *var_params);
};

#endif // VARIATIONALMC_H
