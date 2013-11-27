#ifndef FUNCTION_H
#define FUNCTION_H

/* Class for function implementation.
 *
 * This class is made to implement functions as
 * objects (functor) to be used by polymorphism.
 *
 * Functions are implemented as a sub-class of
 * this class, where the methods are overloaded
 * in order to get the right functionality of
 * the function.
 *
 * Functions can depend on multiple arguments and
 * parameters.
 *
 * Virtual methods without implementations must be
 * implemented when new fuction object is created.
 *
 * Methods:
 *          - contructor  : Takes dimensionality of the function.
 *          - set_params  : Used to set function parameters.
 *          - operator()  : Used to call the function. Must be
 *                          implemented, no standard implementation.
 *          - derivative  : Used to find the first derivative of the
 *                          function.
 *          - derivative2 : Used to find the second derivative of
 *                          the function.
 * */
class Function
{
private:
    double h;
    double *arg_p, *arg_m;

protected:
    int dimension;
    double *params;

public:
    Function(int dimension);

    virtual double set_params(double *params);
    virtual double operator()(double *args){}

    virtual void derivative(double *args, double *diff1);
    virtual double derivative2(double *args);

    virtual ~Function(){}
};

#endif // FUNCTION_H
