/*
Description:

This is a C++ function that takes in two 4x4 matrices X1 and X2 and 
returns the translation error between the two matrices. It makes use 
of the Eigen C++ library, specifically the Vector3d and MatrixXd classes.

The function first extracts the translation vectors p1 and p2 from 
the matrices X1 and X2 respectively, by using the block function to 
select the 3x1 submatrices located at (0,3) within each matrix. 
The difference between these two vectors is then computed using 
the norm function, which calculates the Euclidean norm of the 
resulting 3x1 vector.

The function then returns the computed translation error as a double value.
*/

#include <iostream>
#include <eigen3/Eigen/Dense>

double tranError(Eigen::MatrixXd X1, Eigen::MatrixXd X2) {
    Eigen::Vector3d p1 = X1.block<3,1>(0,3);
    Eigen::Vector3d p2 = X2.block<3,1>(0,3);

    return (p1 - p2).norm();
}

// TEST

/*
int main() {
    Eigen::MatrixXd X1(4, 4);
    X1 << 1, 0, 0, 2,
          0, 1, 0, 3,
          0, 0, 1, 4,
          0, 0, 0, 1;

    Eigen::MatrixXd X2(4, 4);
    X2 << 1, 0, 0, 5,
          0, 1, 0, 6,
          0, 0, 1, 7,
          0, 0, 0, 1;

    double error = tranError(X1, X2);
    std::cout << "Translation error: " << error << std::endl;

    return 0;
}
*/