#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "Function.h"
#include "Integral.h"
#include "UnixTime.h"
#include "problem_definitions.h"
#include "output_functions.h"

using namespace std;

int main(int argc, char *argv[])
{
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
    Integral *integrate;

    a = 1;
    integrand.set_params(&a);

    if( strcmp(method, "GaussLegendre") == 0 ){
        integrate = new GaussLegendre(6);
    }
    if( strcmp(method, "GaussHermite") == 0 ){
        integrate = new GaussHermite(6);
    }
    if( strcmp(method, "MonteCarloBF") == 0 ){
        integrate = new MonteCarloBF(6);
    }
    if( strcmp(method, "MonteCarloIS") == 0 ){
        integrate = new MonteCarloIS(6);
    }

    integral = (*integrate)(lower, upper, N, &integrand);


    if( string(method).find(string("MonteCarlo"))!=string(method).npos )
    {
        variance = (*integrate).get_variance();
        cout << variance << endl;
    }

    cout << integral << endl;

    //test_integrators();
}
