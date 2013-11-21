#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "Function.h"

class Integral
{
protected:
    double integral;
    double *args;

    int dimension;

    Function (*func);

public:
    Integral(){}

    virtual double operator()(double lower, double upper,
                              int n_points, Function f){}

    virtual ~Integral(){}
};


class GaussQuad: public Integral
{
    private:
        void dimension_loops(int N, double *args,
                             int ind, int *indices);

    protected:
        double *x, *w;

        virtual void get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points){}
        virtual double new_term(double *args, int ind, int *indices){}

    public:
        GaussQuad(int dimension);
        double operator()(double lower, double upper,
                          int n_points, Function *f);
};


class GaussLegendre: public GaussQuad
{
    private:
        void get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points);
        double new_term(double *args, int ind, int *indices);

    public:
        GaussLegendre(int dimension);
};

class GaussHermite: public GaussQuad
{
    private:
        void get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points);
        double new_term(double *args, int ind, int *indices);

    public:
        GaussHermite(int dimension);
};

class MonteCarloBF: public Integral
{
    private:
        double integral, variance, term, jacobidet;
        double *args;

    public:
        MonteCarloBF(int dimension);
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
        MonteCarloIS(int dimension);
        double operator()(int n_points, Function *f);
        double get_variance();
};

#endif // INTEGRAL_H
