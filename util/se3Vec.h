/*
DESCRIPTION:

The program defines a function named "se3Vec" that takes a 4x4 matrix as input, 
which represents an element of se(3), and returns a 6x1 vector that represents 
the same element in the Lie algebra of se(3). The function first checks the 
dimension of the input matrix to determine whether it is a skew-symmetric 
matrix or a vector. If the input matrix is skew-symmetric, the function computes 
and returns the corresponding vector. Otherwise, if the input is a vector, 
the function computes and returns the corresponding skew-symmetric matrix. 
The output vector or matrix is represented by an Eigen object of type 
"Eigen::Matrix<double, 6, 1>".
*/

#include <iostream>
#include <eigen3/Eigen/Core>

Eigen::Matrix<double, 6, 1> se3Vec(const Eigen::Matrix4d& X)
{
    Eigen::Matrix<double, 6, 1> g;

    if (X.cols() == 4) // If input is skew-sym change to vector
    {
        g << -X(1,2), X(0,2), -X(0,1), X(0,3), X(1,3), X(2,3);
    } 
    else // If input is vector change to skew-sym
    {       
        g << 0, -X(2), X(1), X(3), 
            X(2), 0, -X(0), X(4), 
            -X(1), X(0), 0, X(5), 
            0, 0, 0, 0; 
    }

    return g; 
} 

/* int main()
{
    Eigen::Matrix<double, 4, 4> X;
    X << 0, -1, 0, 1,
         1, 0, 0, 2,
         0, 0, 1, 3,
         0, 0, 0, 1;
    
    Eigen::Matrix<double, 6, 1> g = se3Vec(X);
    std::cout << "g = \n" << g << "\n";

    return 0;
} */