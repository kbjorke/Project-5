#ifndef FUNCTION_H
#define FUNCTION_H

class Function
{
private:
    double h;

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
