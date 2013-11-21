#include <cmath>
#include <cstdlib>

#include "Integral.h"
#include "Function.h"
#include "lib.h"
#include "hermite.h"
#include "gaussiandeviate.h"


GaussQuad::GaussQuad(int dimension)
{
    this->dimension = dimension;
}



void GaussQuad::dimension_loops(int N, double *args, int ind, int *indices)
{
    int i;

    if( ind == dimension )
    {
        integral += new_term(args, ind, indices);
    }
    else
    {
        for( i = 0; i < N; i++ )
        {
            args[ind] = x[i];
            indices[ind] = i;
            dimension_loops(N, args, ind+1, indices);
        }
    }
}

double GaussQuad::operator()(double lower, double upper,
                  int n_points, Function *f)
{
    args = new double[dimension];
    x = new double[n_points];
    w = new double[n_points];

    int *indices = new int[dimension];

    func = f;

    get_weigths(lower, upper, x, w, n_points);

    integral = 0;

    dimension_loops(n_points, args, 0, indices);

    return integral;
}



GaussLegendre::GaussLegendre(int dimension) : GaussQuad(dimension){}

void GaussLegendre::get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points)
{
    gauleg(lower, upper, x, w, n_points);
}

double GaussLegendre::new_term(double *args, int ind, int *indices)
{
    static int i;
    static double term;

    term = (*func)(args);
    for(i = 0; i < ind; i++ )
    {
        term *= w[indices[i]];
    }
    return term;
}


GaussHermite::GaussHermite(int dimension)
{
    this->dimension = dimension;
}

void GaussHermite::dimension_loops(int N, double *args, int ind, int *indices)
{
    int i;
    double mu;

    if( ind == dimension )
    {

        mu = 0;

        term = (*func)(args);
        for(i = 0; i < ind; i++ )
        {
            term *= w[indices[i]];
            mu += args[i]*args[i];
        }
        // Check for posibility to get exponential part outside
        // dimension_loop, with intention of making GaussQuad
        // Class.
        integral += term*exp(mu);
    }
    else
    {
        for( i = 0; i < N; i++ )
        {
            args[ind] = x[i];
            indices[ind] = i;
            dimension_loops(N, args, ind+1, indices);
        }
    }
}

double GaussHermite::operator()(int n_points, Function *f)
{
    args = new double[dimension];
    x = new double[n_points];
    w = new double[n_points];

    int *indices = new int[dimension];

    func = f;

    gausshermite(x, w, n_points);

    integral = 0;

    dimension_loops(n_points, args, 0, indices);

    return integral;
}

MonteCarloBF::MonteCarloBF(int dimension)
{
    this->dimension = dimension;
}

double MonteCarloBF::operator()(double lower, double upper,
                                int n_points, Function *f)
{
    int i, dim;
    double random_num;

    args = new double[dimension];
    func = f;

    integral = 0;
    variance = 0;

    jacobidet = 1;
    for( dim = 0; dim < dimension; dim++ )
    {
        jacobidet *= (upper - lower);
    }

    srand(time(NULL));
    for( i = 0; i < n_points; i++ )
    {
        for( dim = 0; dim < dimension; dim++ )
        {
            random_num = (double)rand() /  RAND_MAX;
            args[dim] = lower + random_num*(upper - lower);
        }

        term = (*func)(args);
        integral += term;
        variance += term*term;
   }

    integral = integral/((double) n_points);
    variance = variance/((double) n_points) - integral*integral;

    integral = jacobidet*integral;
    variance = jacobidet*sqrt(variance/((double) n_points));

    return integral;
}

double MonteCarloBF::get_variance()
{
    return variance;
}


MonteCarloIS::MonteCarloIS(int dimension)
{
    this->dimension = dimension;
}

double MonteCarloIS::operator()(int n_points, Function *f)
{
    int i, dim;
    double random_num, mu, sqrt2;

    long int idum = -1;

    args = new double[dimension];
    func = f;

    integral = 0;
    variance = 0;

    jacobidet = pow(acos(-1),dimension/2); //1;
    /*
    for( dim = 0; dim < dimension; dim++ )
    {
        jacobidet *= (upper - lower);
    }
    */

    sqrt2 = 1/sqrt(2);

    srand(time(NULL));
    for( i = 0; i < n_points; i++ )
    {
        mu = 0;
        for( dim = 0; dim < dimension; dim++ )
        {
            random_num = gaussian_deviate(&idum)*sqrt2;
            args[dim] = random_num;
            mu += args[dim]*args[dim];
        }

        term = (*func)(args)*exp(mu);
        integral += term;
        variance += term*term;
   }

    integral = integral/((double) n_points);
    variance = variance/((double) n_points) - integral*integral;

    integral = jacobidet*integral;
    variance = jacobidet*sqrt(variance/((double) n_points));

    return integral;
}

double MonteCarloIS::get_variance()
{
    return variance;
}
