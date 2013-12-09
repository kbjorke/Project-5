#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "mpi.h"
#include "globals.h"
#include "Function.h"
#include "Integral.h"
#include "Hamiltonian.h"
#include "variational_MC.h"
#include "MetropolisQM.h"
#include "UnixTime.h"
#include "problem_definitions.h"
#include "output_functions.h"


using namespace std;

int main(int argc, char *argv[])
{
    double integral, upper, lower, variance, a, energy, delta;
    double stdev, acceptance_rate, delta_t;
    double *R_init, *params, *var_params;
    int i, N, ind_var, num_params, num_var;
    char *method, *trial_psi;
    string *id_params, trial_lable;
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

    //if( solve_integral && strcmp(method, "Metropolis") != 0 )
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

    /*
    else if( solve_integral && strcmp(method, "Metropolis") == 0 )
    {
        Function *psi_trial;

        psi_trial = new Psi_T1(6);

        params = new double[1];
        *params = 1;
        (*psi_trial).set_params(params);

        Observable *potential_energy;
        potential_energy = new PotentialEnergy();

        MetropolisQM metroQM(6);

        delta = 0.7;

        metroQM.initialize(psi_trial, delta);

        if( my_rank == 0 )
        {
            t0 = getUnixTime();
        }

        energy = metroQM(potential_energy, N);

        if( my_rank == 0 )
        {
            t1 = getUnixTime();

            stdev = metroQM.get_std();
            acceptance_rate = metroQM.get_acceptance_rate();
            delta_t = t1 - t0;

            cout << energy*pow(acos(-1),3) << endl;
            cout << stdev << endl;
            cout << acceptance_rate << endl;
            cout << delta_t << endl;

        }
    }
    */

    else if( integrators_test )
    {
        test_integrators();
    }

    else if( VMC )
    {
        Function *psi_trial;

        if( strcmp(trial_psi, "T1") == 0 ){
            psi_trial = new Psi_T1(6);

            trial_lable = "Psi_T1";

            num_params = 1;
            num_var = 10001;
            params = new double[num_params];

            *params = 1;

            id_params =  new string[num_params];
            id_params[0] = "alpha";

            ind_var = 0;

            var_params =  new double[num_var];
            for( i = 0; i < num_var; i++ )
            {
                var_params[i] = 0.5 + i*0.0001;
            }
        }
        if( strcmp(trial_psi, "T2") == 0 ){
            psi_trial = new Psi_T2(6);

            trial_lable = "Psi_T2";

            num_params = 2;
            num_var = 1001;

            params = new double[num_params];

            params[0] = 1.0;

            id_params =  new string[num_params];
            id_params[0] = "alpha";
            id_params[1] = "beta";

            ind_var = 1;

            var_params =  new double[num_var];
            for( i = 0; i < num_var; i++ )
            {
                var_params[i] = 0.2 + i*0.0005;
            }
        }

        Function *potential;

        if( nointeract )
        {
            potential = new V_HO(6);

            trial_lable.append(" no-interaction");
        }
        else
        {
            potential = new V_TOT(6);

            trial_lable.append(" interaction");
        }


        Observable *hamiltonian;

        hamiltonian = new Hamiltonian();

        (*hamiltonian).set_potential(potential);

        double best_var;

        energy = variational_MC(psi_trial, hamiltonian,
                                params, num_params, var_params,
                                ind_var, num_var, &best_var,
                                N, id_params, true, trial_lable);

        if( my_rank == 0 )
        {
        cout << energy << endl;
        }
    }

    MPI_Finalize();
}
