#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <ctime>

#include <iostream>

#include "mpi.h"
#include "globals.h"
#include "output_functions.h"
#include "problem_definitions.h"
#include "Integral.h"
#include "MetropolisQM.h"
#include "UnixTime.h"

using namespace std;

void test_integrators()
{
    int i, j;

    string methods[4];

    methods[0] = "GaussLegendre";
    methods[1] = "GaussHermite";
    methods[2] = "MonteCarloBF";
    methods[3] = "MonteCarloIS";

    output_method(methods[3]);
    /*
    for( i = 2; i < 4; i++ )
    {
        output_method(methods[i]);
    }
    */
}

void output_method(string method)
{
    int i, j;
    long int N;
    double a, lower, upper, integral, variance, t0, t1, delta_t;

    Integral *integrate;
    Integrand integrand(6);
    a = 1;
    integrand.set_params(&a);

    lower = -4;
    upper = 4;

    // Generate filename for output file

    fstream outfile;

    if( my_rank == 0 )
    {
        string output_file;
        ostringstream oss;

        oss << "output_" << method  << "_np" << numprocs << ".txt";

        output_file = oss.str();

        outfile.open(output_file.c_str(), ios::out);
    }

    if( method.find("Gauss")!=method.npos )
    {
        if( method.find("GaussLegendre")!=method.npos )
        {
            integrate = new GaussLegendre(6);
        }
        if( method.find("GaussHermite")!=method.npos )
        {
            integrate = new GaussHermite(6);
        }

        if( my_rank == 0)
        {
            outfile << "Method: " << method << '\t' <<
                       "Processors: " << numprocs << endl;
            outfile << "N" << '\t' << "N*dim" << '\t' <<
                       "Integral-value" << '\t' <<
                       "Time" << endl;

            cout << method << endl;
            outfile.precision(10);
        }

        for( N = 6; N <= 60; N += 6 )
        {
            if( my_rank == 0 )
            {
                cout << N << endl;
                t0 = getUnixTime();
            }

            integral = (*integrate)(lower, upper, N, &integrand);

            if( my_rank == 0 )
            {
                t1 = getUnixTime();

                delta_t = t1-t0;

                outfile << setw(2) << setfill(' ') << N << '\t' <<
                           setw(11) << setfill(' ') << pow(N,6) << '\t' <<
                           integral << '\t' <<
                           delta_t << endl;
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }

    }
    if( method.find("MonteCarlo")!=method.npos )
    {
        if( method.find("MonteCarloBF")!=method.npos )
        {
            integrate = new MonteCarloBF(6);
        }
        if( method.find("MonteCarloIS")!=method.npos )
        {
            integrate = new MonteCarloIS(6);
        }

        if( my_rank == 0 )
        {
            outfile << "Method: " << method << '\t' <<
                       "Processors: " << numprocs << endl;
            outfile << '\t' <<
                       "N" << '\t' << "Integral-value" << '\t' <<
                       "Variance" << '\t' << "Time" << endl;

            cout << method << endl;
            outfile.precision(10);
        }

        for( N = 100; N <= 1e10; N *= 10 )
        {
            if( my_rank == 0 )
            {
                cout << N << endl;
                t0 = getUnixTime();
            }

            integral = (*integrate)(lower, upper, N, &integrand);

            if( my_rank == 0 )
            {
                t1 = getUnixTime();

                delta_t = t1-t0;

                variance = (*integrate).get_variance();

                outfile << setw(11) << setfill(' ') << N << '\t' <<
                           integral << '\t' <<
                           variance << '\t' <<
                           time << endl;
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    if( my_rank == 0 )
    {
        outfile.close();
    }
}

void output_VMC_data_header(fstream *outfile, int N, int ind_var,
                            double *params, int num_params,
                            string *id_params)
{
    static int i;

    string output_file;
    ostringstream oss;

    char datetime[80];

    time_t rawtime;
    tm* timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(datetime, 80, "%F-%H%M", timeinfo);

    oss << "output_VMC_" << "_np" << numprocs << "_" <<
           datetime << ".txt";

    output_file = oss.str();

    (*outfile).open(output_file.c_str(), ios::out);

    (*outfile) << "QM Variational Monte Carlo: " << endl;
    (*outfile) << "Processors: " << numprocs << '\t' <<
                  "N: " << N << endl;
    (*outfile) << "Variational parameter: " << id_params[ind_var] << endl;
    (*outfile) << "Non-variational parameters: " << endl;

    if( num_params != 1 )
    {
        for( i = 0; i < num_params; i++ )
        {
            if( i != ind_var )
            {
                (*outfile) << id_params[i] << " = " <<
                              params[i] << " ," <<'\t';
            }
        }
    }
    else
    {
        (*outfile) << "No non-variational parameters.";
    }

    (*outfile) << endl << setw(20) << setfill('-') << " " << endl;

    (*outfile) << id_params[ind_var] << '\t' <<
                  "Energy" << '\t' <<
                  "Standard deviance" << '\t' <<
                  "Acceptance rate" << '\t' <<
                  "Time" << endl;

    (*outfile).precision(10);
}

void output_VMC_data(fstream *outfile, double *params,
                     int ind_var, double energy, double stdev,
                     double acceptance_rate, double delta_t)
{
    (*outfile) << params[ind_var] << '\t' <<
                  energy << '\t' <<
                  stdev << '\t' <<
                  acceptance_rate << '\t' <<
                  delta_t << endl;
}

void output_VMC_data_end(fstream *outfile, double best_energy,
                         string *id_params, int ind_best,
                         double best_var)
{
        (*outfile) << setw(20) << setfill('-') << " " << endl;
        (*outfile) << "Best energy: " << best_energy << '\t' <<
                      " with " << id_params[ind_best] <<
                      " = " << best_var << endl;

        (*outfile).close();
}
