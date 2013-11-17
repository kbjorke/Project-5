#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Function.h"
#include "Integral.h"

using namespace std;

class Integrand: public Function
{
    private:
    double alpha, mu, rho;
    double dist[3];

    public:
    Integrand(double alpha)
    {
        this->alpha = alpha;
    }


    double operator()(double *r)
    {
        int i;

        mu = 0;
        for( i = 0; i < 3; i++ )
        {
            mu += r[0+i]*r[0+i] + r[3+i]*r[3+i];

            dist[i] = r[3+i] - r[i];
        }

        rho = sqrt(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);

        return exp(-alpha*alpha*mu)/rho;
    }
};

class Func: public Function
{
    private:

    public:
    Func(){}


    double operator()(double *r)
    {
        //return exp(-(*r))/(*r);
        return 1/(2 + (*r)*(*r));
    }
};

int main(int argc, char *argv[])
{
    double integral;
    Func func;

    GaussLegendre integrate(1);

    integral = integrate(0,3,10,func);

    cout << integral << endl;
}
