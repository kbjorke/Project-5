#include <cmath>
#include <cstdlib>

#include "Integral.h"
#include "lib.h"
#include "hermite.h"


GaussLegendre::GaussLegendre(double dimension)
{
    this->dimension = dimension;
}

void GaussLegendre::dimension_loops(int N, double *param, int ind, int *indices)
{
    int i;

    if( ind == dimension )
    {
        term = (*func)(param);
        for(i = 0; i < ind; i++ )
        {
            term *= w[indices[i]];
        }
        integral += term;
    }
    else
    {
        for( i = 0; i < N; i++ )
        {
            param[ind] = x[i];
            indices[ind] = i;
            dimension_loops(N, param, ind+1, indices);
        }
    }
}

double GaussLegendre::operator()(double lower, double upper,
                  int n_points, Function *f)
{
    param = new double[dimension];
    x = new double[n_points];
    w = new double[n_points];

    int *indices = new int[dimension];

    func = f;

    gauleg(lower, upper, x, w, n_points);

    integral = 0;

    dimension_loops(n_points, param, 0, indices);

    return integral;
}

GaussHermite::GaussHermite(double dimension)
{
    this->dimension = dimension;
}

void GaussHermite::dimension_loops(int N, double *param, int ind, int *indices)
{
    int i;
    double mu;

    if( ind == dimension )
    {

        mu = 0;

        term = (*func)(param);
        for(i = 0; i < ind; i++ )
        {
            term *= w[indices[i]];
            mu += param[i]*param[i];
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
            param[ind] = x[i];
            indices[ind] = i;
            dimension_loops(N, param, ind+1, indices);
        }
    }
}

double GaussHermite::operator()(int n_points, Function *f)
{
    param = new double[dimension];
    x = new double[n_points];
    w = new double[n_points];

    int *indices = new int[dimension];

    func = f;

    gausshermite(x, w, n_points);

    integral = 0;

    dimension_loops(n_points, param, 0, indices);

    return integral;
}

MonteCarloBF::MonteCarloBF(double dimension)
{
    this->dimension = dimension;
}

double MonteCarloBF::operator()(double lower, double upper,
                                int n_points, Function *f)
{
    int i, dim;
    double random_num;

    param = new double[dimension];
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
            param[dim] = lower + random_num*(upper - lower);
        }

        term = (*func)(param);
        integral += term;
        variance += term*term;
    }

    integral = jacobidet*integral/((double) n_points);
    variance = jacobidet*variance/((double) n_points);

    return integral;
}

double MonteCarloBF::get_variance()
{
    return variance;
}
