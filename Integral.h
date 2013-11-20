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
        double *args, *x, *w;

        void dimension_loops(int N, double *args, int ind, int *indices);

    public:

        GaussLegendre(double dimension);
        double operator()(double lower, double upper,
                          int n_points, Function *f);
};

class GaussHermite: public Integral
{
    private:
        double integral, term;
        double *args, *x, *w;

        void dimension_loops(int N, double *args, int ind, int *indices);

    public:

        GaussHermite(double dimension);
        double operator()(int n_points, Function *f);
};

class MonteCarloBF: public Integral
{
    private:
        double integral, variance, term, jacobidet;
        double *args;

    public:
        MonteCarloBF(double dimension);
        double operator()(double lower, double upper,
                          int n_points, Function *f);
        double get_variance();
};

class MonteCarloIS: public Integral
{
    private:
        double integral, variance, term, jacobidet;
        double *args;

    public:
        MonteCarloIS(double dimension);
        double operator()(double lower, double upper,
                          int n_points, Function *f);
        double get_variance();
};

#endif // INTEGRAL_H
