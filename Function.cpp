#include "Function.h"

#include <iostream>

using namespace std;

Function::Function(int dimension)
{
    this->dimension = dimension;

    h = 0.0001;
}

void Function::derivative(double *args, double *diff1)
{
    int i;
    static double arg, arg_p, arg_m;
    diff1 = new double[dimension];

    for( i = 0; i < dimension; i++ )
    {
        arg = *args;
        arg_p = arg + h;
        arg_m = arg - h;

        diff1[i] = ((*this)(&arg_p) - (*this)(&arg_m))/(2*h);
        cout << diff1[i] << endl;
    }
}


double Function::derivative2(double *args)
{
    int i;
    static double arg, arg_p, arg_m;
    double diff2 = 0;

    for( i = 0; i < dimension; i++ )
    {
        arg = *args;
        arg_p = arg + h;
        arg_m = arg - h;

        diff2 += ((*this)(&arg_p) - 2*(*this)(&arg) + (*this)(&arg_m))/(h*h);
        cout << (*this)(&arg_m) << endl;
    }
    return diff2;
}
