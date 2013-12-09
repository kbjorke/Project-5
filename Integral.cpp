#include <cmath>
#include <cstdlib>

#include "mpi.h"
#include "globals.h"
#include "Integral.h"
#include "Function.h"
#include "lib.h"
#include "hermite.h"
#include "gaussiandeviate.h"
#include "UnixTime.h"

/* Constructor for Integral superclass.
 *
 * Takes the dimensionality of the integral as input and
 * allocates memory for arguments array.
 *
 * Input:
 *          dimension : Dimensionality of the integral,
 *                      i.e. the number of variables to
 *                      be given as argument to the function.
 * */
Integral::Integral(int dimension)
{
    this->dimension = dimension;
    args = new double[dimension];
}

/* Constructor for Gaussian quadrature class.
 *
 * Takes the dimensionality of the integral as passes
 * on the the Integral superclass, and allocates
 * memory to the indices array, to be used in dimension loops.
 * */
GaussQuad::GaussQuad(int dimension) : Integral(dimension)
{
    indices = new int[dimension];
}

/* Method dimension_loops.
 *
 * Used to preform a nested loop over all the dimensions
 * in order to set arguments for the function evaluation.
 *
 * Input:
 *                N : Number of points the integral will be
 *                    evaluated per dimension.
 *            *args : Pointer or array to be filled with
 *                    arguments from the nested loops.
 *              ind : The current index of the argument, or
 *                    the current level of the loop.
 *         *indices : Pointer or array containing the index
 *                    of the current loop.
 * */
void GaussQuad::dimension_loops(int N, double *args, int ind, int *indices)
{
    int i;

    // Enter when the loop has reached its innermost level:
    if( ind == dimension )
    {
        integral += new_term(args, indices);
    }
    else
    {
        for( i = 0; i < N; i++ )
        {
            // Sets argument-points for the current dimension.
            args[ind] = x[i];
            indices[ind] = i;
            dimension_loops(N, args, ind+1, indices);
        }
    }
}

/* Method operator().
 *
 * Used as a functor to solve the integral for a given function,
 * integration points and integration limits.
 *
 * Input:
 *           lower : Lower integration limit.
 *           upper : Upper integration limit.
 *        n_points : Number of points for the evaluation
 *                   of the integral.
 *              *f : Pointer to a functor object which
 *                   represents the integrand.
 * */
double GaussQuad::operator()(double lower, double upper,
                  int n_points, Function *f)
{
    int i;
    double final_integral;

    x = new double[n_points];
    w = new double[n_points];

    get_weigths(lower, upper, x, w, n_points);

    func = f;

    integral = 0;
    final_integral = 0;

    for( i = my_rank; i < n_points; i += numprocs )
    {
        args[0] = x[i];
        indices[0] = i;
        dimension_loops(n_points, args, 1, indices);
    }

    MPI_Reduce(&integral, &final_integral, 1, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);


    if( my_rank == 0 )
    {
        return final_integral;
    }
    else
    {
        return 0;
    }

    delete[] x;
    delete[] w;
}


GaussLegendre::GaussLegendre(int dimension) : GaussQuad(dimension){}

/* Method get_weigths.
 *
 * Used the function gauleg from the library lib.cpp to get the
 * weigthts and integration points for Gauss-Legendre integration.
 *
 * Input:
 *           lower : Lower integration limit.
 *           upper : Upper integration limit.
 *              *x : Pointer or array to be filled with the
 *                   integration points for gaussian-legendre
 *                   integration.
 *              *w : Pointer or array to be filled with the
 *                   integration weigths for gaussian-legendre
 *                   integration.
 *        n_points : Number of points for the evaluation
 *                   of the integral.
 * */
void GaussLegendre::get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points)
{
    gauleg(lower, upper, x, w, n_points);
}

/* Method new_term.
 *
 * Used to calculate new term for the integral.
 *
 * Input:
 *           *args : Pointer or array which contains arguments
 *                   for where the integrand will be evaluated.
 *        *indices : Pointer or array containing indexes for
 *                   current loop which are used t get the
 *                   right weigth.
 * */
double GaussLegendre::new_term(double *args, int *indices)
{
    static int i;
    static double term;

    term = (*func)(args);
    for(i = 0; i < dimension; i++ )
    {
        term *= w[indices[i]];
    }
    return term;
}


GaussHermite::GaussHermite(int dimension) : GaussQuad(dimension){}

/* Method get_weigths.
 *
 * Used the function gausshermite to get the weigthts and integration
 * points for Gauss-Legendre integration.
 *
 * Input:
 *           lower : Lower integration limit.
 *           upper : Upper integration limit.
 *              *x : Pointer or array to be filled with the
 *                   integration points for gaussian-legendre
 *                   integration.
 *              *w : Pointer or array to be filled with the
 *                   integration weigths for gaussian-legendre
 *                   integration.
 *        n_points : Number of points for the evaluation
 *                   of the integral.
 * */
void GaussHermite::get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points)
{
    gausshermite(x, w, n_points);
}

/* Method new_term.
 *
 * Used to calculate new term for the integral. Multiplies the
 * term with exp(mu) wich is the inverse of the weigthfunction
 * for Gauss-Hermite integration.
 *
 * Input:
 *           *args : Pointer or array which contains arguments
 *                   for where the integrand will be evaluated.
 *        *indices : Pointer or array containing indexes for
 *                   current loop which are used t get the
 *                   right weigth.
 * */
double GaussHermite::new_term(double *args, int *indices)
{
    static int i;
    static double term, mu;

    mu = 0;

    term = (*func)(args);
    for(i = 0; i < dimension; i++ )
    {
        term *= w[indices[i]];
        mu += args[i]*args[i];
    }

    return term*exp(mu);
}



/* Constructor for Gaussian quadrature class.
 *
 * Takes the dimensionality of the integral as passes
 * on the the Integral superclass.
 * */
MonteCarlo::MonteCarlo(int dimension) : Integral(dimension){}

/* Method operator().
 *
 * Used as a functor to solve the integral for a given function,
 * integration points and integration limits.
 *
 * Input:
 *           lower : Lower integration limit.
 *           upper : Upper integration limit.
 *        n_points : Number of points for the evaluation
 *                   of the integral.
 *              *f : Pointer to a functor object which
 *                   represents the integrand.
 * */
double MonteCarlo::operator()(double lower, double upper,
                                int n_points, Function *f)
{
    static long int time_int;
    int i;
    double local_integral, local_variance;

    func = f;

    integral = 0;
    variance = 0;
    local_integral = 0;
    local_variance = 0;

    constant = constant_term(lower, upper);

    //In order to make sure that the seed can change between quick calls
    time_int = (long int) (getUnixTime()*10000 + my_rank);

    set_seed(time_int);

    for( i = my_rank; i < n_points; i += numprocs )
    {
        term = new_term(lower, upper);
        local_integral += term;
        local_variance += term*term;
    }

    MPI_Reduce(&local_integral, &integral, 1, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);
    MPI_Reduce(&local_variance, &variance, 1, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);

    if( my_rank == 0 )
    {
        integral = integral/((double) n_points);
        variance = variance/((double) n_points) - integral*integral;

        integral = constant*integral;
        variance = constant*sqrt(variance/((double) n_points));

        return integral;
    }
    else
    {
        return 0;
    }
}

/* Method get_variance.
 *
 * Return the variance for the previously calcualted integral.
 * */
double MonteCarlo::get_variance()
{
    return variance;
}

MonteCarloBF::MonteCarloBF(int dimension) : MonteCarlo(dimension){}

/* Method constant_term.
 *
 * Returns constant term to be used in integration, here the
 * Jacobideterminant because of change of variables.
 *
 * Input:
 *           lower : Lower integration limit.
 *           upper : Upper integration limit.
 * */
double MonteCarloBF::constant_term(double upper, double lower)
{
    static int dim;
    static double jacobidet;

    jacobidet = 1;
    for( dim = 0; dim < dimension; dim++ )
    {
        jacobidet *= (upper - lower);
    }
    return jacobidet;
}

/* Method new_term.
 *
 * Return new term for the integral, evaluated by a
 * uniformly distributed random position.
 *
 * Input:
 *           lower : Lower integration limit.
 *           upper : Upper integration limit.
 * */
double MonteCarloBF::new_term(double lower, double upper)
{
    static int dim;
    static double random_num;

    for( dim = 0; dim < dimension; dim++ )
    {
        random_num = ran2(&idum);
        args[dim] = lower + random_num*(upper - lower);
    }

    return (*func)(args);
}

/* Method set_seed.
 *
 * Sets seed for random number generator.
 *
 * Input:
 *          - seed : Seed to be set.
 * */
void MonteCarloBF::set_seed(long int seed)
{
    static long int seed_;

    // To make sure that the seed is not too large for the
    // ran2() random number function:
    seed_ = seed - ((long int) seed/(int) 1e7)*(int) 1e7;
    idum = -seed_;
}


MonteCarloIS::MonteCarloIS(int dimension) : MonteCarlo(dimension){}

/* Method constant_term.
 *
 * Returns constant term to be used in integration, here
 * from a change of variables.
 *
 * Input:
 *           lower : Lower integration limit.
 *           upper : Upper integration limit.
 * */
double MonteCarloIS::constant_term(double upper, double lower)
{
    return pow(acos(-1),dimension/2);
}

/* Method new_term.
 *
 * Return new term for the integral, evaluated by a
 * gaussian distributed random position. Multilpies
 * the function evaluation with exp(mu), which is the
 * inverse of the probability distribution for
 * Importance Sampling with a gaussian PDF.
 *
 * Input:
 *           lower : Lower integration limit.
 *           upper : Upper integration limit.
 *
 * NOTE: This is a to general implementation, should be rewritten
 *       so that is can take a arbritary PDF.
 * */
double MonteCarloIS::new_term(double lower, double upper)
{
    static int dim;
    static double random_num, mu;

    mu = 0;
    for( dim = 0; dim < dimension; dim++ )
    {
        random_num = gaussian_deviate(&idum)/sqrt(2);
        args[dim] = random_num;
        mu += args[dim]*args[dim];
    }

    return (*func)(args)*exp(mu);
}

/* Method set_seed.
 *
 * Sets seed for random number generator.
 *
 * Input:
 *          - seed : Seed to be set.
 * */
void MonteCarloIS::set_seed(long int seed)
{
    static long int seed_;

    // To make sure that the seed is not too large for the
    // ran2() random number function:
    seed_ = seed - ((long int) seed/(int) 1e7)*(int) 1e7;
    idum = -seed_;
}
