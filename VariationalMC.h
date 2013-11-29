#ifndef VARIATIONALMC_H
#define VARIATIONALMC_H

#include "Function.h"
#include "Hamiltonian.h"

class VariationalMC
{
private:
    int dimension, accepts, N, thermalization;
    double delta, mean_energy, std_energy;
    double energy, variance;
    double *var_params, *R, *R_new;
    Hamiltonian *H;
    Function *psi;

    int numprocs, my_rank;

    double get_localenergy(double *R);
    void get_newpos(double *R, double *R_new);
    void accept_test(double *R, double *R_new);
    double relative_probability(double *R, double *R_new);
    void set_seed(double seed);

public:
    VariationalMC(int dimension);
    void initialize(Hamiltonian *H, double delta, double *R_init);
    double operator()(Function *trial_psi, double *var_params, int n_points);
    double get_std();
    double get_acceptance_rate();

    ~VariationalMC(){}
};

#endif // VARIATIONALMC_H
