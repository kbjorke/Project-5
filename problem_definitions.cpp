#include <cmath>

#include "problem_definitions.h"

#include <iostream>

using namespace std;


Integrand::Integrand(int dimension) : Function(dimension){}

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


Psi_T1::Psi_T1(int dimension) : Function(dimension){}

double Psi_T1::set_params(double *params)
{
    alpha = *params;
}

double Psi_T1::operator()(double *r)
{
    int i;

    mu = 0;
    for( i = 0; i < 3; i++ )
    {
        mu += r[0+i]*r[0+i] + r[3+i]*r[3+i];
    }

    return exp(-alpha*alpha*mu*0.5);
}


Psi_T2::Psi_T2(int dimension) : Function(dimension){}

double Psi_T2::set_params(double *params)
{
    alpha = params[0];
    beta = params[1];
}

double Psi_T2::operator()(double *r)
{
    int i;

    mu = 0;
    for( i = 0; i < 3; i++ )
    {
        mu += r[0+i]*r[0+i] + r[3+i]*r[3+i];

        dist[i] = r[3+i] - r[i];
    }
    rho = sqrt(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);

    exponent = -alpha*alpha*mu*0.5 + rho*0.5/(1 + beta*rho);

    return exp(exponent);
}

V_TOT::V_TOT(int dimension) : Function(dimension){}

double V_TOT::operator()(double *r)
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
        return 0.5*mu + 1/rho;
    }
}


V_HO::V_HO(int dimension) : Function(dimension){}

double V_HO::operator()(double *r)
{
    int i;

    mu = 0;
    for( i = 0; i < 3; i++ )
    {
        mu += r[0+i]*r[0+i] + r[3+i]*r[3+i];
    }

    return 0.5*mu;
}
