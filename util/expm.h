/*
DESCRIPTION:

The given code implements the matrix exponential function using the 
Taylor series expansion method. The function takes an input matrix A 
of type Eigen::MatrixXd and returns its matrix exponential, exp(A), 
of the same type.

The function first initializes expA to the identity matrix of the same 
size as A. It then initializes a variable factorial to 1 and a variable
i to 1. It also initializes An to A.

The function then enters a do-while loop that calculates the terms of 
the Taylor series expansion of exp(A) using the formula exp(A) = I + A/1! + A^2/2! + A^3/3! + ... 
and updates expA and An at each iteration. The variable i is used to 
calculate the factorials needed for each term.

The loop exits when the norm of the current term of the Taylor series 
expansion, An/factorial, is smaller than a certain threshold, which is 
set to 1e-15 in this case.

Finally, the function returns the matrix expA, which should be an 
approximation of exp(A) calculated using the Taylor series expansion method.
*/

#include <cmath>
#include <eigen3/Eigen/Dense>

Eigen::MatrixXd expm(const Eigen::MatrixXd& A) 
{
    Eigen::MatrixXd expA = Eigen::MatrixXd::Identity(A.rows(), A.cols());
    double factorial = 1;
    int i = 1;
    Eigen::MatrixXd An = A;
    
    do
    {
        factorial *= i;
        expA += An / factorial;
        An = An * A;
        ++i;
    } while ((An / factorial).lpNorm<Eigen::Infinity>() > 1e-15);
  
    return expA;
}