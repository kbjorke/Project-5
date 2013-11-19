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
    Integrand(double alpha)
    {
        this->alpha = alpha;
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
    Func(){}

    double operator()(double *r)
    {
        //double x;

        //x = (*r);

        //return x*x*exp(-x*x - 2*x);
        return (*r)*(*r)*exp(-(*r)*(*r) - 2*(*r));
    }
};

int main(int argc, char *argv[])
{
    double integral, upper, lower;
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

    /*
    Integrand integrand(1);

    if( strcmp(method, "GaussLegendre") == 0 ){
        GaussLegendre integrate(6);
        integral = integrate(lower, upper, N, &integrand);
    }
    if( strcmp(method, "GaussHermite") == 0 ){
        GaussHermite integrate(6);
        integral = integrate(N, &integrand);
    }
    */

    Func func;
    GaussHermite integrate(1);
    cout << integrate(N, &func);

    // Want to have integral evaluation here, problems with integral
    // variable, not declared in scope.
    // integral = integrate(lower, upper, N, &integrand);

    cout << integral << endl;
}
