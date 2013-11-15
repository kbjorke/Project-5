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


    double operator()(double **r)
    {
        mu = r[0][0]*r[0][0] + \
                r[0][1]*r[0][1] + \
                r[0][2]*r[0][2] + \
                r[1][0]*r[1][0] + \
                r[1][1]*r[1][1] + \
                r[1][2]*r[1][2];

        dist[0] = r[1][0] - r[0][0];
        dist[1] = r[1][1] - r[0][1];
        dist[2] = r[1][2] - r[0][2];

        rho = sqrt(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]);

        return exp(-alpha*alpha*mu)/rho;
    }
};



int main(int argc, char *argv[])
{
}
