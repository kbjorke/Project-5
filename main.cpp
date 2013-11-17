#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Function.h"

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

void inner_loops(int maxloop, int maxvalue, int ind, int *val)
{
    int i;

    if( ind == maxloop )
    {
        for(i = 0; i < ind; i++ )
        {
            cout << val[i];
        }
        cout << endl;
    }
    else
    {
        for( i = 0; i < maxvalue; i++ )
        {
            val[ind] = i;
            inner_loops(maxloop, maxvalue, ind+1, val);
        }
    }
}

void multi_loops(int maxloop, int maxvalue)
{
    int i;
    int *val = new int[maxloop];

    inner_loops(maxloop, maxvalue, 0, val);
}


int main(int argc, char *argv[])
{
    multi_loops(5,2);
}
