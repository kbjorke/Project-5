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

    for( i = 0; i < 4; i++ )
    {
        output_method(methods[i]);
    }
    /*
    string output_file;
    ostringstream oss;

    // Generate filename for output file
    oss << "data_" << method  << "_n" <<
           n << "_T" <<
           setw(3) << setfill('0') << int(T*100) << "_alpha" <<
           setw(5) << setfill('0') << int(round(alpha*10000)) <<
           ".txt" << '\0' << endl;

    output_file = oss.str();

    fstream myfile;
    myfile.open(output_file.c_str(), ios::out);

    // Write out first line in output file, general data about solution
    myfile << "T: " << T << '\t' <<
              "d: " << d << '\t' <<
              "dt: " << delta_t << '\t' <<
              "dx: " << delta_x << '\t' <<
              "n: " << n << '\t' <<
              "m: " << m << '\t' <<
              "alpha: " << alpha << '\t' <<
              "method: " << method << '\t' << endl;

    // Writes out solution
    myfile.precision(10);
    myfile << scientific;
    for( j = 0; j < m; j++ ){
        for( i = 0; i < n; i++ )
        {
            myfile << u[i][j] << "   ";
        }
        myfile << endl;
    }

    myfile.close();
    */
}

void output_method(string method)
{
    int i, j, N;
    double a, lower, upper, integral, t0, t1, time;

    string output_file;
    ostringstream oss;

    Integral *integrate;
    Integrand integrand(6);
    a = 1;
    integrand.set_params(&a);


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

        outfile << "Method: " << method << '\t' << endl;

        for( N = 6; N <= 20; N += 2 )
        {
            t0 = getUnixTime();
            //integral = integrate(lower, upper, N, &integrand);
            t1 = getUnixTime();

            time = t1-t0;

            outfile << N << '\t' <<
                       pow(N,6) << '\t' <<
                       integral << '\t' <<
                       time << endl;
        }

    }
    if( method.find("MonteCarlo")!=method.npos )
    {
        if( method.find("MonteCarloBF")!=method.npos )
        {
            MonteCarloBF integrate(6);
        }
        if( method.find("MonteCarloIS")!=method.npos )
        {
            MonteCarloIS integrate(6);
        }

        outfile << "Method: " << method << '\t' << endl;
    }

    /*
    fstream myfile;
    myfile.open(output_file.c_str(), ios::out);

    // Write out first line in output file, general data about solution
    myfile << "T: " << T << '\t' <<
              "d: " << d << '\t' <<
              "dt: " << delta_t << '\t' <<
              "dx: " << delta_x << '\t' <<
              "n: " << n << '\t' <<
              "m: " << m << '\t' <<
              "alpha: " << alpha << '\t' <<
              "method: " << method << '\t' << endl;

    // Writes out solution
    myfile.precision(10);
    myfile << scientific;
    for( j = 0; j < m; j++ ){
        for( i = 0; i < n; i++ )
        {
            myfile << u[i][j] << "   ";
        }
        myfile << endl;
    }
    */
    outfile.close();
}


