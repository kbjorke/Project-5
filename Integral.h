#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "Function.h"

class Integral
{
protected:
    double *params;
    int dimension;

    Function (*func);

public:
    Integral(){}

    virtual double operator()(double lower, double upper,
                              int n_points, Function f){}

    virtual ~Integral(){}
};



class GaussLegendre: public Integral
{
    private:

        double integral, term;
        double *param, *x, *w;

        void dimension_loops(int N, double *param, int ind, int *indices);

    public:

        GaussLegendre(double dimension);
        double operator()(double lower, double upper,
                          int n_points, Function *f);
};

class GaussHermite: public Integral
{
    private:
        double integral, term;
        double *param, *x, *w;

        void dimension_loops(int N, double *param, int ind, int *indices);

    public:

        GaussHermite(double dimension);
        double operator()(int n_points, Function *f);
};

#endif // INTEGRAL_H
