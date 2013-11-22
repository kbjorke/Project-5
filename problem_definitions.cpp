#include <cmath>

#include "problem_definitions.h"


Integrand::Integrand(int dimension) : Function(dimension)
{
}

double Integrand::set_params(double *params)
{
    alpha = *params;
}

double Integrand::operator()(double *r)
{
    int i;

    mu = 0;
    for( i = 0; i < 3; i++ )
    {
        mu += r[0+i]*r[0+i] + r[3+i]*r[3+i];

        dist[i] = r[3+i] - r[i];
    }
    rho = sqrt(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);

    if( rho == 0 )
    {
        return 0;
    }
    else
    {
        return exp(-alpha*alpha*mu)/rho;
    }
}


