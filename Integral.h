#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "Function.h"
#include <iostream>

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
                              int n_points, Function *f){}
    virtual double get_variance(){}

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



class MonteCarlo: public Integral
{
    protected:
        double variance, term, constant;

        virtual double constant_term(double upper, double lower){}
        virtual double new_term(double upper, double lower){}
        virtual void set_seed(double seed){}

    public:
        MonteCarlo(int dimension);
        double operator()(double lower, double upper,
                          int n_points, Function *f);
        double get_variance();
};


class MonteCarloBF: public MonteCarlo
{
    public:
        MonteCarloBF(int dimension);
        double constant_term(double upper, double lower);
        double new_term(double lower, double upper);
        void set_seed(double seed);
};


class MonteCarloIS: public MonteCarlo
{
    private:
        long int idum;
    public:
        MonteCarloIS(int dimension);
        double constant_term(double upper, double lower);
        double new_term(double upper, double lower);
        void set_seed(double seed);
};

#endif // INTEGRAL_H
