#ifndef METROPOLISQM_H
#define METROPOLISQM_H

#include "Function.h"
#include "Observable.h"

class MetropolisQM
{
private:
    int dimension, accepts, N, total_accepts;
    double delta, mean_energy, std_energy;
    double energy, variance;
    double *var_params, *R, *R_new;
    Observable *O;
    Function *psi;

    double get_localenergy(double *R);
    void get_newpos(double *R, double *R_new);
    void accept_test(double *R, double *R_new);
    double relative_probability(double *R, double *R_new);
    void set_seed(double seed);

public:
    MetropolisQM(int dimension);
    void initialize(Function *psi, double delta, double *R_init);
    double operator()(Observable *O, int n_points);
    double get_std();
    double get_acceptance_rate();

    ~MetropolisQM(){}
};

#endif // METROPOLISQM_H
