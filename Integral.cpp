#include <cmath>
#include <cstdlib>

#include "globals.h"
#include "mpi.h"
#include "Integral.h"
#include "Function.h"
#include "lib.h"
#include "hermite.h"
#include "gaussiandeviate.h"
#include "UnixTime.h"

Integral::Integral(int dimension)
{
    this->dimension = dimension;
    args = new double[dimension];
}

GaussQuad::GaussQuad(int dimension) : Integral(dimension)
{
    indices = new int[dimension];
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

    delete[] x;
    delete[] w;
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


GaussHermite::GaussHermite(int dimension) : GaussQuad(dimension){}

void GaussHermite::get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points)
{
    gausshermite(x, w, n_points);
}

double GaussHermite::new_term(double *args, int ind, int *indices)
{
    static int i;
    static double term, mu;

    mu = 0;

    term = (*func)(args);
    for(i = 0; i < ind; i++ )
    {
        term *= w[indices[i]];
        mu += args[i]*args[i];
    }

    return term*exp(mu);
}



MonteCarlo::MonteCarlo(int dimension) : Integral(dimension){}

double MonteCarlo::operator()(double lower, double upper,
                                int n_points, Function *f)
{
    int i;
    double local_integral, local_variance;

    func = f;

    integral = 0;
    variance = 0;
    local_integral = 0;
    local_variance = 0;

    constant = constant_term(lower, upper);

    set_seed((int) getUnixTime()*100+my_rank);

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
}

double MonteCarlo::get_variance()
{
    return variance;
}

MonteCarloBF::MonteCarloBF(int dimension) : MonteCarlo(dimension){}

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

double MonteCarloBF::new_term(double lower, double upper)
{
    static int dim;
    static double random_num;

    for( dim = 0; dim < dimension; dim++ )
    {
        random_num = ran2(&idum);
        //random_num = (double)rand() /  RAND_MAX;
        args[dim] = lower + random_num*(upper - lower);
    }

    return (*func)(args);
}

void MonteCarloBF::set_seed(int seed)
{
    idum = -seed;
    //srand(seed);
}


MonteCarloIS::MonteCarloIS(int dimension) : MonteCarlo(dimension){}

double MonteCarloIS::constant_term(double upper, double lower)
{
    return pow(acos(-1),dimension/2);
}

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

void MonteCarloIS::set_seed(int seed)
{
    idum = -seed;
}
