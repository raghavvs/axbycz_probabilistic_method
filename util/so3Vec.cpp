/*
DESCRIPTION:

The above program defines a function named so3_vec that takes a 3D vector 
as input and returns a 3x3 matrix. If the input is a vector, the function 
converts it to a skew-symmetric matrix, and if the input is a skew-symmetric 
matrix, the function returns it directly. The resulting 3x3 matrix represents 
an element of the special orthogonal group SO(3), which is used in 3D 
rotation calculations.
*/

#include <iostream>
#include <eigen3/Eigen/Core>

Eigen::MatrixXd so3Vec(const Eigen::MatrixXd& X)
{
    if (X.cols() == 3) // If input is skew-sym change to vector
    {
        Eigen::Matrix<double, 3, 1> g;
        g << -X(1,2), X(0,2), -X(0,1);
        return g;
    }
    else // If input is vector change to skew-sym
    {
        Eigen::Matrix3d g;
        g << 0, -X(2), X(1),
            X(2), 0, -X(0),
            -X(1), X(0), 0;
        return g;
    }
}

int main()
{
    Eigen::Matrix<double, 3, 3> M;
    M << 0, -3, 2,
        3, 0, -1,
        -2, 1, 0;

    Eigen::Matrix<double, 3, 1> vec = so3Vec(M);
    std::cout << "Mat2Vec = \n" << vec << "\n";

    Eigen::Matrix<double, 3, 1> V;
    V << 1, 2, 3;

    Eigen::Matrix3d mat = so3Vec(V);
    std::cout << "Vec2Mat = \n" << mat << "\n";

    return 0;
}