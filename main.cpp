#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "Function.h"
#include "Integral.h"

using namespace std;

class Integrand: public Function
{
    private:
    double alpha, mu, rho;
    double dist[3];

    public:
    Integrand(int dimension) : Function(dimension)
    {
    }

    double set_params(double *params)
    {
        alpha = *params;
    }

    double operator()(double *r)
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
};

class Func: public Function
{
    private:

    public:
    Func(int dimension) : Function(dimension)
    {
    }

    double operator()(double *r)
    {
        return 2*(*r)*(*r) + 4*(*r) + 6;
    }
};


int main(int argc, char *argv[])
{
    /*
    double integral, upper, lower, variance, a;
    int i, N;
    char *method;

    // Loop over commandline arguments to find parameters and options:
    for( i = 0; i < argc-1; i++ ){
        if( strcmp(argv[i], "-method") == 0 ){
            method = argv[i+1];
        }
        if( strcmp(argv[i], "-lower") == 0 ){
            lower = atof(argv[i+1]);
        }
        if( strcmp(argv[i], "-upper") == 0 ){
            upper = atof(argv[i+1]);
        }
        if( strcmp(argv[i], "-N") == 0 ){
            N = atoi(argv[i+1]);
        }
        if( strcmp(argv[i], "-N") == 0 ){
            N = atoi(argv[i+1]);
        }
    }

    Integrand integrand(6);

    a = 1;
    integrand.set_params(&a);

    if( strcmp(method, "GaussLegendre") == 0 ){
        GaussLegendre integrate(6);
        integral = integrate(lower, upper, N, &integrand);
    }
    if( strcmp(method, "GaussHermite") == 0 ){
        GaussHermite integrate(6);
        integral = integrate(N, &integrand);
    }
    if( strcmp(method, "MonteCarloBF") == 0 ){
        MonteCarloBF integrate(6);
        integral = integrate(lower, upper, N, &integrand);
        variance = integrate.get_variance();
        cout << variance << endl;
    }
    if( strcmp(method, "MonteCarloIS") == 0 ){
        MonteCarloIS integrate(6);
        integral = integrate(N, &integrand);
        variance = integrate.get_variance();
        cout << variance << endl;
    }

    // Want to have integral evaluation here, problems with integral
    // variable, not declared in scope.
    // integral = integrate(lower, upper, N, &integrand);

    cout << integral << endl;
    */

    Func func(1);

    double x = 0;
    double *diff1 = new double[1];

    func.derivative(&x, diff1);

    cout << func(&x) << endl;
    cout << diff1[0] << endl;
    cout << func.derivative2(&x) << endl;
}
