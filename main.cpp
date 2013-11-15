#include <iostream>
#include <cstdlib>

#include "Function.h"

using namespace std;

class func2x: public Function
{
    public:
    func2x(){}

    double operator()(double x)
    {
        return 2*x;
    }
};

class funcxy: public Function
{
    public:
    funcxy(){}

    double x,y;

    double operator()(double *variables)
    {
        x = variables[0];
        y = variables[1];

        cout << x << " " << y << endl;

        return x*y;
    }
};



int main(int argc, char *argv[])
{
    func2x func1;
    funcxy func2;
    double x;
    double *variables;

    x = atof(argv[argc-2]);

    variables = new double[2];
    variables[0] = x;
    variables[1] = atof(argv[argc-1]);

    cout << func1(x) << endl;
    cout << func2(variables) << endl;
}
