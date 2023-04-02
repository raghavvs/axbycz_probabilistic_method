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

double metric(const std::vector<Eigen::Matrix4d>& A,
              const std::vector<Eigen::Matrix4d>& B,
              const std::vector<Eigen::Matrix4d>& C,
              const Eigen::Matrix4d& X,
              const Eigen::Matrix4d& Y,
              const Eigen::Matrix4d& Z) {
    double diff = 0.0;
    int N = 0;

    for (size_t i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].cols() / 4; ++j) { // divide by 4 to avoid out-of-bounds block access
            Eigen::Matrix4d lhs = A[i].block(0, j * 4, 4, 4) * X *
                                  B[i].block(0, j * 4, 4, 4);
            Eigen::Matrix4d rhs = Y * C[i].block(0, j * 4, 4, 4) * Z;
            diff += (lhs - rhs).norm();
            N++;
        }
    }

    diff /= static_cast<double>(N);
    return diff;
}

int main() {
    // Set the number of A, B, and C matrices
    int num_matrices = 3;

    // Create random 4x4 matrices for A, B, C, X, Y, and Z
    std::vector<Eigen::Matrix4d> A(num_matrices);
    std::vector<Eigen::Matrix4d> B(num_matrices);
    std::vector<Eigen::Matrix4d> C(num_matrices);
    Eigen::Matrix4d X = Eigen::Matrix4d::Random();
    Eigen::Matrix4d Y = Eigen::Matrix4d::Random();
    Eigen::Matrix4d Z = Eigen::Matrix4d::Random();

    for (int i = 0; i < num_matrices; ++i) {
        A[i] = Eigen::Matrix4d::Random();
        B[i] = Eigen::Matrix4d::Random();
        C[i] = Eigen::Matrix4d::Random();
    }

    // Compute the metric value
    double metric_value = metric(A, B, C, X, Y, Z);

    // Output the result
    std::cout << "Metric value: " << metric_value << std::endl;

    return 0;
}