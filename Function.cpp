#include "Function.h"

#include <iostream>

using namespace std;

Function::Function(int dimension)
{
    int i;

    this->dimension = dimension;

    arg_m = new double[dimension];
    arg_p = new double[dimension];

    h = 0.0001;
}

double Function::set_params(double *params)
{
    this->params = params;
}

void Function::derivative(double *args, double *diff1)
{
    static int i, j;

    for( i = 0; i < dimension; i++ )
    {
        for( j = 0; j < dimension; j++ )
        {
            if( i == j )
            {
                arg_p[i] = args[i] + h;
                arg_m[i] = args[i] - h;
            }
            else
            {
                arg_p[i] = args[i];
                arg_m[i] = args[i];
            }
        }

        diff1[i] = ((*this)(arg_p) - (*this)(arg_m))/(2*h);
        cout << diff1[i] << endl;
    }
}


double Function::derivative2(double *args)
{
    static int i, j;
    static double diff2 = 0;

    for( i = 0; i < dimension; i++ )
    {
        for( j = 0; j < dimension; j++ )
        {
            if( i == j )
            {
                arg_p[i] = args[i] + h;
                arg_m[i] = args[i] - h;
            }
            else
            {
                arg_p[j] = args[j];
                arg_m[j] = args[j];
            }
        }

        //cout << (*this)(arg_p) - 2*(*this)(args) + (*this)(arg_m) << endl;
        diff2 += ((*this)(arg_p) - 2*(*this)(args) + (*this)(arg_m))/(h*h);
    }
    return diff2;
}
