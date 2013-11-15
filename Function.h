#ifndef FUNCTION_H
#define FUNCTION_H

class Function
{
public:
    Function();

    virtual double operator()(double x){}
    virtual double operator()(double *x){}
    virtual double operator()(double **x){}

    virtual ~Function(){}
};

#endif // FUNCTION_H
