/*
DESCRIPTION:

The program computes the Frobenius norm of a difference between two matrix 
expressions, each composed of matrices from multiple 3D arrays. The difference 
is computed as the product of three matrices (A, X, and B) minus the product 
of three other matrices (Y, C, and Z), where A, B, C, X, Y, and Z are matrix 
arrays. The program iterates over all matrices in the third dimension of each 
matrix array and computes the norm of the difference between the two matrix 
expressions for each matrix. The sum of the norms is divided by the total 
number of matrices to compute the average norm.
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>

double metric(const std::vector<Eigen::Matrix3d>& A,
              const std::vector<Eigen::Matrix3d>& B,
              const std::vector<Eigen::Matrix3d>& C,
              const Eigen::Matrix3d& X,
              const Eigen::Matrix3d& Y,
              const Eigen::Matrix3d& Z) {
    double diff = 0.0;
    int N = 0;
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].rows(); ++j) {
            for (int k = 0; k < A[i].cols(); ++k) {
                diff += std::abs((A[i](j, k) * X(k, j) * B[i](j, k) - Y(j, k) * C[i](k, j) * Z(j, k)));
                ++N;
            }
        }
    }
    return diff / N;
}

/* int main()
{
    std::cout << "Build successful" << std::endl;
    return 0;
}  */