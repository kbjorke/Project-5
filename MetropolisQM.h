#ifndef METROPOLISQM_H
#define METROPOLISQM_H

#include "Function.h"
#include "Observable.h"

/* Class for Metropolis algorithm for Quantum Mechanical integration
 *
 * Used to preform the Metropolis algorithm on a wavefunction and
 * an observable to find the expectation value of the observable
 * in the wavefunction.
 *
 *
 * Methods:
 *    - contructor           : Takes dimensionality of the integral.
 *    - initialize           : Set the wavefunction *psi, and the step
 *                             length, delta, of the random walk.
 *    - operator()           : Used to solve for the expectation-value
 *                             of a given observable O.
 *    - get_std              : Returns the standard deviation of the
 *                             latest calculated expectation value.
 *    - get_acceptance_rate  : Return the rate of how often the new
 *                             steps was accepted.
 *    - set_seed             : Set the seed for the random number generator.
 *    - set_thermalization   : Set the thermalizaion (how many steps we want
 *                             the random walker to take before we begin to
 *                             sample.
 *
 *    - get_localeigenvalue  : Get the local eigenvalue of the observable,
 *                             to be stored as a sample.
 *    - get_newpos           : Gets a new step from the previous step and the
 *                             step length.
 *    - accept_test          : Test if the new step is accepted or not.
 *    - relative_probability : Used to calculate the relative probability
 *                             probability between the new position and
 *                             the previous.
 *
 *
 * */
class MetropolisQM
{
private:
    int dimension, accepts, N, total_accepts, thermalization;
    double delta, mean_energy, std_energy;
    double energy, variance;
    double *var_params, *R, *R_new;
    Observable *O;
    Function *psi;

    long int idum;

    double get_localeigenvalue(double *R);
    void get_newpos(double *R, double *R_new);
    void accept_test(double *R, double *R_new);
    double relative_probability(double *R, double *R_new);

public:
    MetropolisQM(int dimension);
    void initialize(Function *psi, double delta);
    double operator()(Observable *O, int n_points);
    double get_std();
    double get_acceptance_rate();
    void set_seed(long int seed);
    void set_thermalization(int thermalization);

    ~MetropolisQM(){}
};

#endif // METROPOLISQM_H
