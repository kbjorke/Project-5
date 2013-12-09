#include "Function.h"

using namespace std;

/* Constructor for Function class.
 *
 * Takes the dimensionality of the function as input and
 * sets some nececary variables like the argument arrays
 * and h-variable used for the derivation methods.
 *
 * Input:
 *          dimension : Dimensionality of the function,
 *                      i.e. the number of variables to
 *                      be given as argument to the function.
 * */
Function::Function(int dimension)
{
    int i;

    this->dimension = dimension;

    arg_m = new double[dimension];
    arg_p = new double[dimension];

    h = 0.0001; // To be used by numerical derivative.
}

/* Method set_params.
 *
 * Used to set the parameters used by the function.
 * This is a general implementation, which can be used
 * in a function implementation, but one can also make
 * a specific implementation since this is a virtual
 * method.
 *
 * Input:
 *          *params : Pointer or array containing the
 *                    parameters to be set.
 * */
double Function::set_params(double *params)
{
    this->params = params;
}

/* Method derivative.
 *
 * General implementation of the first-derivative/gradient of the
 * function, based on a 2 point formula with a steplength h which
 * is set in the contructor.
 *
 * In the case where a closed for expression for the first derivative
 * exist, this could be implemented by overwriting this virtual method.
 *
 * Input:
 *          *args  : Pointer or array containing the arguments
 *                   (position) where the derivative will be evaluated.
 *
 *          *diff1 : Pointer or array where the derivative/gradient
 *                   is to be stored.
 *
 * NOTE: This method could be improved or split up to get a implementation
 *      	 of the derivative of the function in respect to a given variable.
 * */
void Function::derivative(double *args, double *diff1)
{
    static int i, j;

    // Loop to set the forward and backward points.
    for( i = 0; i < dimension; i++ )
    {
        for( j = 0; j < dimension; j++ )
        {
            // Change for the given variable.
            if( i == j )
            {
                arg_p[i] = args[i] + h;
                arg_m[i] = args[i] - h;
            }
            // Hold the other variables as they where.
            else
            {
                arg_p[i] = args[i];
                arg_m[i] = args[i];
            }
        }

        // Find the derivative in respect to the given variable.
        diff1[i] = ((*this)(arg_p) - (*this)(arg_m))/(2*h);
    }
}

/* Method derivative2.
 *
 * General implementation of the second-derivative/laplacian of the
 * function, based on a 3 point formula with a steplength h which
 * is set in the contructor.
 *
 * In the case where a closed for expression for the second derivative
 * exist, this could be implemented by overwriting this virtual method.
 *
 * Input:
 *          *args  : Pointer or array containing the arguments
 *                   (position) where the derivative will be evaluated.
 * */
double Function::derivative2(double *args)
{
    static int i, j;
    static double diff2;

    diff2 = 0;
    for( i = 0; i < dimension; i++ )
    {
        for( j = 0; j < dimension; j++ )
        {
            // Change for the given variable.
            if( i == j )
            {
                arg_p[i] = args[i] + h;
                arg_m[i] = args[i] - h;
            }
            // Hold the other variables as they where.
            else
            {
                arg_p[j] = args[j];
                arg_m[j] = args[j];
            }
        }

        // Find the derivative in respect to the given variable.
        diff2 += ((*this)(arg_p) - 2*(*this)(args) + (*this)(arg_m))/(h*h);
    }
    return diff2;
}
