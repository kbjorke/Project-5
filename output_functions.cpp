#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>

#include <iostream>

#include "output_functions.h"
#include "problem_definitions.h"
#include "Integral.h"
#include "UnixTime.h"

using namespace std;

/* Function that writes solution to file.
 *
 * Input:
 *
 *			- n: amount of points to describe position
 *
 *			- m: amount of point to describe time
 *
 *			- T: final time of solution
 *
 *			- d: largest positions
 *
 *			- delta_t: time step
 *
 *			- delta_x: spatial step
 *
 *			- alpha: alpha value for solution, given by delta_t/delta_x^2
 *
 *			- *method: name of methods used for solution
 *
 * 			- **u: n x m - matrix containing solution
 *
 *
 * Output:
 *
 *			- Filename: data_<method>_n<n>_T<T>_alpha<alpha>.txt
 *
 * */
void test_integrators()
{
    int i, j;

    string methods[4];

    methods[0] = "GaussLegendre";
    methods[1] = "GaussHermite";
    methods[2] = "MonteCarloBF";
    methods[3] = "MonteCarloIS";

    //output_method(methods[2]);
    for( i = 0; i < 4; i++ )
    {
        output_method(methods[i]);
    }
}

void output_method(string method)
{
    int i, j;
    long int N;
    double a, lower, upper, integral, variance, t0, t1, time;

    string output_file;
    ostringstream oss;

    Integral *integrate;
    Integrand integrand(6);
    a = 1;
    integrand.set_params(&a);

    lower = -4;
    upper = 4;

    // Generate filename for output file
    oss << "output_" << method  << ".txt";

    output_file = oss.str();

    fstream outfile;
    outfile.open(output_file.c_str(), ios::out);

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

        outfile << "Method: " << method << endl;
        outfile << "N" << '\t' << "N*dim" << '\t' << "Integral-value" << '\t' <<
                   "Time" << endl;

        cout << method << endl;
        outfile.precision(10);
        for( N = 6; N <= 60; N += 6 )
        {
            cout << N << endl;
            t0 = getUnixTime();
            integral = (*integrate)(lower, upper, N, &integrand);
            t1 = getUnixTime();

            time = t1-t0;

            outfile << setw(2) << setfill(' ') << N << '\t' <<
                       setw(11) << setfill(' ') << pow(N,6) << '\t' <<
                       integral << '\t' <<
                       time << endl;
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

        outfile << "Method: " << method << endl;
        outfile << '\t' <<
                   "N" << '\t' << "Integral-value" << '\t' <<
                   "Variance" << '\t' << "Time" << endl;

        cout << method << endl;
        outfile.precision(10);
        for( N = 100; N <= 1e10; N *= 10 )
        {
            cout << N << endl;
            t0 = getUnixTime();
            integral = (*integrate)(lower, upper, N, &integrand);
            t1 = getUnixTime();

            time = t1-t0;

            variance = (*integrate).get_variance();

            outfile << setw(11) << setfill(' ') << N << '\t' <<
                       integral << '\t' <<
                       variance << '\t' <<
                       time << endl;
        }
    }

    outfile.close();
}


