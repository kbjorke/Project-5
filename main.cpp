#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "globals.h"
#include "mpi.h"
#include "Function.h"
#include "Integral.h"
#include "Hamiltonian.h"
#include "MetropolisQM.h"
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
    char *method, *trial_psi;
    bool solve_integral = false;
    bool integrators_test = false;
    bool VMC = false;
    bool nointeract = false;

    double t0, t1, time;

    // Loop over commandline arguments to find parameters and options:
    for( i = 0; i < argc; i++ ){
        if( strcmp(argv[i], "-method") == 0 ){
            method = argv[i+1];
        }
        if( strcmp(argv[i], "-trial") == 0 ){
            trial_psi = argv[i+1];
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
        if( strcmp(argv[i], "-integral") == 0 ){
            solve_integral = true;
        }
        if( strcmp(argv[i], "-int_test") == 0 ){
            integrators_test = true;
        }
        if( strcmp(argv[i], "-VMC") == 0 ){
            VMC = true;
        }
        if( strcmp(argv[i], "-nointeract") == 0 ){
            nointeract = true;
        }
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    if( solve_integral )
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

        if( my_rank == 0 )
        {
            t0 = getUnixTime();
        }
        integral = (*integrate)(lower, upper, N, &integrand);
        if( my_rank == 0 )
        {
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

    else if( integrators_test )
    {
        test_integrators();
    }

    else if( VMC )
    {
        Function *psi_trial;

        if( strcmp(trial_psi, "T1") == 0 ){
            psi_trial = new Psi_T1(6);
        }
        if( strcmp(trial_psi, "T2") == 0 ){
            psi_trial = new Psi_T1(6);
        }

        Function *potential;

        if( nointeract )
        {
            potential = new V_HO(6);
        }
        else
        {
            potential = new V_TOT(6);
        }


        Observable *hamiltonian;

        hamiltonian = new Hamiltonian();

        (*hamiltonian).set_potential(potential);

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

        alpha = 1;

        (*psi_trial).set_params(&alpha);

        MetropolisQM metroQM(6);

        metroQM.initialize(psi_trial, delta, R_init);


        if( my_rank == 0 )
        {
            t0 = getUnixTime();
        }

        energy = metroQM(hamiltonian, N);

        if( my_rank == 0 )
        {
            t1 = getUnixTime();

            std = metroQM.get_std();
            acceptance_rate = metroQM.get_acceptance_rate();
            time = t1 - t0;

            cout << energy << endl;
            cout << std << endl;
            cout << acceptance_rate << endl;
            cout << time << endl;
        }
    }

    MPI_Finalize();
}
