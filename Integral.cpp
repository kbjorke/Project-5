#include <cmath>
#include <cstdlib>
#include <mpi.h>

#include "Integral.h"
#include "Function.h"
#include "lib.h"
#include "hermite.h"
#include "gaussiandeviate.h"
#include "UnixTime.h"
//#include "globals.h"


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
    int i, N;
    double h, final_integral;

    int numprocs, my_rank;
    int argc;
    char **argv;

    args = new double[dimension];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if( my_rank == 0 )
    {
        N = n_points/numprocs;
        h = (upper-lower)/numprocs;
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&h, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    integral = my_rank;

    x = new double[N];
    w = new double[N];

    int *indices = new int[dimension];

    func = f;

    get_weigths(lower + h*my_rank, upper - h*(numprocs-my_rank),
                x, w, N);

    integral = 0;

    dimension_loops(N, args, 0, indices);

    if( my_rank == 0 )
    {
        MPI_Status status;
        final_integral = integral;
        for( i = 0; i < numprocs; i++ )
        {
            MPI_Recv(&integral, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 500,
                     MPI_COMM_WORLD, &status);
            final_integral += integral;
        }
        cout << final_integral << endl;

        return final_integral;
    }
    else
    {
        MPI_Send(&integral, 1, MPI_DOUBLE, 0, 500, MPI_COMM_WORLD);
    }

    MPI_Finalize();
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



MonteCarlo::MonteCarlo(int dimension)
{
    this->dimension = dimension;
}

double MonteCarlo::operator()(double lower, double upper,
                                int n_points, Function *f)
{
    int i, dim;

    args = new double[dimension];
    func = f;

    integral = 0;
    variance = 0;

    constant = constant_term(lower, upper);

    for( i = 0; i < n_points; i++ )
    {
        term = new_term(lower, upper);
        integral += term;
        variance += term*term;
    }

    integral = integral/((double) n_points);
    variance = variance/((double) n_points) - integral*integral;

    integral = constant*integral;
    variance = constant*sqrt(variance/((double) n_points));

    return integral;
}

double MonteCarlo::get_variance()
{
    return variance;
}

MonteCarloBF::MonteCarloBF(int dimension) : MonteCarlo(dimension)
{
    srand(getUnixTime()*1000);
}

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
        random_num = (double)rand() /  RAND_MAX;
        args[dim] = lower + random_num*(upper - lower);
    }

    return (*func)(args);
}


MonteCarloIS::MonteCarloIS(int dimension) : MonteCarlo(dimension)
{
    idum = -(getUnixTime()*100);
}

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
