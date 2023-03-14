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

Eigen::MatrixXd se3Vec(const Eigen::MatrixXd& X)
{
    if (X.cols() == 4) // If input is skew-sym change to vector
    {
        Eigen::Matrix<double, 6, 1> g;
        g << -X(1,2), X(0,2), -X(0,1), X(0,3), X(1,3), X(2,3);
        return g;
    } 
    else // If input is vector change to skew-sym
    {
        Eigen::Matrix4d g;
        g << 0, -X(2), X(1), X(3), 
            X(2), 0, -X(0), X(4), 
            -X(1), X(0), 0, X(5), 
            0, 0, 0, 0;
        return g;
    }
} 

 int main()
{
    Eigen::Matrix<double, 4, 4> M;
    M << 0, -3, 2, 4,
         3, 0, -1, 5,
         -2, 1, 0, 6,
         0, 0, 0, 0;

    Eigen::Matrix<double, 6, 1> vec = se3Vec(M);
    std::cout << "Mat2Vec = \n" << vec << "\n";

    Eigen::Matrix<double, 6, 1> V;
    V << 1, 2, 3, 4, 5, 6;

    Eigen::Matrix4d mat = se3Vec(V);
    std::cout << "Vec2Mat = \n" << mat << "\n";

    return 0;
}