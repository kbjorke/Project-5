#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "Function.h"

/* Superclass for Integral implementations.
 *
 * This class is made to implement integrals as
 * subclasses.
 *
 * The methods and variable in this class are things
 * which all integrals have in common.
 *
 * Most of the methods must be implemented in the next
 * classlevel in a way which characterize that specific
 * type of integrators.
 *
 *
 * Methods:
 *          - contructor    : Takes dimensionality of the integral.
 *          - operator()    : Used to solve the integral for a given
 *                            integrand, number of points and integration
 *                            limits. Must be implemented, no standard
 *                            implementation.
 *          - get_variance  : Used to return the variance of the
 *                            previously solved integra.
 *
 * */
class Integral
{
protected:
    double integral;
    double *args;

    int dimension;

    Function (*func);

public:
    Integral(int dimension);

    virtual double operator()(double lower, double upper,
                              int n_points, Function *f){}
    virtual double get_variance(){}

    virtual ~Integral(){}
};



/* Class for Gaussian Quadrature integral implementation.
 *
 * This class is made to be the superclass of integrals of
 * the type Gaussian Quadrature.
 *
 *
 * Methods:
 *          - contructor     : Takes dimensionality of the integral.
 *          - operator()     : Used to solve the integral. Implemented
 *                             here, and is the same for all Gaussian
 *                             quarature integrals.
 *
 *          - get_weigths    : Virtual method that gets the weigth points
 *                             for the integration, must be implemented
 *                             for each specific GQ implementation.
 *          - new_term       : Virtual method to get the new term for the
 *                             integration. Must be implemented.
 *
 *          - dimension_loops: Used to preform loops over the number of
 *                             dimensions the integral will be solved for.
 *
 * */
class GaussQuad: public Integral
{
    private:
        void dimension_loops(int N, double *args,
                             int ind, int *indices);

    protected:
        int *indices;
        double *x, *w;

        virtual void get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points){}
        virtual double new_term(double *args, int *indices){}

    public:
        GaussQuad(int dimension);
        double operator()(double lower, double upper,
                          int n_points, Function *f);
};


/* Subclass with the implementation of Gauss-Legendre quadrature.
 *
 * Methods:
 *          - contructor     : Takes dimensionality of the integral.
 *
 *          - get_weigths    : Method that gets the gauss-legendre
 *                             weigth points for the integration.
 *          - new_term       : Method to get the new term for the
 *                             integration.
 *
 * */
class GaussLegendre: public GaussQuad
{
    private:
        void get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points);
        double new_term(double *args, int *indices);

    public:
        GaussLegendre(int dimension);
};


/* Subclass with the implementation of Gauss-Hermite quadrature.
 *
 * Methods:
 *          - contructor     : Takes dimensionality of the integral.
 *
 *          - get_weigths    : Method that gets the gauss-hermite
 *                             weigth points for the integration.
 *          - new_term       : Method to get the new term for the
 *                             integration.
 *
 * */
class GaussHermite: public GaussQuad
{
    private:
        void get_weigths(double lower, double upper,
                                 double *x, double *w,
                                 int n_points);
        double new_term(double *args, int *indices);

    public:
        GaussHermite(int dimension);
};



/* Class for Monte Carlo integral implementation.
 *
 * This class is made to be the superclass of integrals of
 * the Monte Carlo type.
 *
 *
 * Methods:
 *          - contructor     : Takes dimensionality of the integral.
 *          - operator()     : Used to solve the integral. Implemented
 *                             here, and is the same for all Monte Carlo
 *                             implementations.
 *          - get_variance   : Get the variance for the solution of a
 *                             Monte Carlo integration.
 *          - set_seed       : Sets the seed for the random number generator
 *                             in use, must be implementet.
 *
 *          - constant_term  : Virtual method that gives a constant value
 *                             to be multiplied to the integral, due to some
 *                             Jacobideterminant, normalization, or
 *                             change of variables, must be implemented
 *                             for each specific MC implementation.
 *          - new_term       : Virtual method to get the new term for the
 *                             integration. Must be implemented.
 *
 * */
class MonteCarlo: public Integral
{
    protected:
        double variance, term, constant;

        virtual double constant_term(double upper, double lower){}
        virtual double new_term(double upper, double lower){}

    public:
        MonteCarlo(int dimension);
        double operator()(double lower, double upper,
                          int n_points, Function *f);
        double get_variance();
        virtual void set_seed(long int seed){}
};


/* Subclass with the implementation of Brute Force Monte Carlo.
 *
 * Methods:
 *          - contructor     : Takes dimensionality of the integral.
 *
 *          - constant_term  : Method that return the constant term
 *                             which comes from the change of
 *                             variables for the integration limits.
 *          - new_term       : Method to get the new term for the
 *                             integration.
 *          - set_seed       : Set the seed for the random numbers.
 *
 * */
class MonteCarloBF: public MonteCarlo
{
    private:
        long int idum;
    public:
        MonteCarloBF(int dimension);
        double constant_term(double upper, double lower);
        double new_term(double lower, double upper);
        void set_seed(long int seed);
};


/* Subclass with the implementation of Importace Sampling Monte Carlo.
 *
 * Methods:
 *          - contructor     : Takes dimensionality of the integral.
 *
 *          - constant_term  : Method that return the constant term
 *                             which comes from the change of
 *                             variables for the random numbers.
 *          - new_term       : Method to get the new term for the
 *                             integration.
 *          - set_seed       : Set the seed for the random numbers.
 *
 * */
class MonteCarloIS: public MonteCarlo
{
    private:
        long int idum;
    public:
        MonteCarloIS(int dimension);
        double constant_term(double upper, double lower);
        double new_term(double upper, double lower);
        void set_seed(long int seed);
};

#endif // INTEGRAL_H
