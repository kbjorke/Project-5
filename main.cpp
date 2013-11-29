#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "Function.h"
#include "Integral.h"
#include "Hamiltonian.h"
#include "VariationalMC.h"
#include "UnixTime.h"
#include "problem_definitions.h"
#include "output_functions.h"

using namespace std;

int main(int argc, char *argv[])
{
    double integral, upper, lower, variance, a;
    double energy, std, acceptance_rate, delta, alpha;
    double *R_init;
    int i, N;
    char *method;
    bool integrators_test = false;
    bool VMC = false;

    double t0, t1, time;

    // Loop over commandline arguments to find parameters and options:
    for( i = 0; i < argc; i++ ){
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
        if( strcmp(argv[i], "-int_test") == 0 ){
         integrators_test = true;
        }
        if( strcmp(argv[i], "-VMC") == 0 ){
         VMC = true;
        }
    }


    if( integrators_test )
    {
        test_integrators();
    }

    else if( VMC )
    {
        V_HO potential(6);
        Hamiltonian H_HO;

        H_HO.set_potential(&potential);

        R_init = new double[6];
        for( i = 0; i < 6; i++ )
        {
            if( i < 3 )
            {
                    R_init[i] = 1.0;
            }
            else
            {
                R_init[i] = -1.0;
            }
        }

        delta = 0.6;

        VariationalMC VMC(6);

        VMC.initialize(&H_HO, delta, R_init);

        Psi_T1 psi_t1(6);

        alpha = 1;

        t0 = getUnixTime();
        energy = VMC(&psi_t1, &alpha, N);
        t1 = getUnixTime();

        std = VMC.get_std();
        acceptance_rate = VMC.get_acceptance_rate();
        time = t1 - t0;

        cout << energy << endl;
        cout << std << endl;
        cout << acceptance_rate << endl;
        cout << time << endl;

    }

    else
    {
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

        t0 = getUnixTime();
        integral = (*integrate)(lower, upper, N, &integrand);
        t1 = getUnixTime();


        if( string(method).find(string("MonteCarlo"))!=string(method).npos )
        {
            variance = (*integrate).get_variance();
            cout << variance << endl;
        }

        time = t1-t0;

        cout << integral << endl;
        cout << time << endl;
    }
}
