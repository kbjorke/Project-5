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
        term = dimension*(*func)(param);
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

    if( ind == dimension )
    {
        term = dimension*(*func)(param);
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
