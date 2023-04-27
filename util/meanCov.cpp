/*
DESCRIPTION:

This program defines several functions for calculating the mean
and covariance of a set of 4x4 matrices. The log function takes
a 4x4 matrix and calculates its matrix logarithm, while the vex
function takes a 3x3 matrix and returns its vector of exteriorization.
The meanCov function takes an array of 3x3 matrices and its size
N and calculates the mean and covariance of the logarithms of those
matrices. It does this by first taking the average of the logarithms
using the expm function, then iteratively refining this average until
convergence using the log and vex functions. Finally, it calculates
the covariance by taking the vector of differences between each logarithm
and the mean logarithm and computing their outer product.

Input:
    X: Matrix dim - 4x4 - pass by reference
    N: Number of A, B, C matrices or data pairs
Output:
    Mean: Matrix dim - 4x4
    Covariance: Matrix dim - 6x6
    No return value to function meanCov - outputs can be obtained from the function parameters
*/

#include <iostream>
#include <cstdlib>
#include <random>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

Eigen::Vector3d vex(const Eigen::Matrix3d &m) {
    Eigen::Vector3d v;
    v << m(2, 1), m(0, 2), m(1, 0);
    return v;
}

void meanCov(const std::vector<Eigen::Matrix4d> &X,
             Eigen::Matrix4d &Mean,
             Eigen::Matrix<double, 6, 6> &Cov) {

    int N = X.size();
    Mean = Eigen::Matrix4d::Identity();
    Cov = Eigen::Matrix<double, 6, 6>::Zero();

    // Initial approximation of Mean
    Eigen::Matrix4d sum_se = Eigen::Matrix4d::Zero();
    for (int i = 0; i < N; i++) {
        sum_se += X[i].log();
    }
    Mean = ((1.0 / N) * sum_se).exp();

    // Iterative process to calculate the true Mean
    Eigen::Matrix4d diff_se = Eigen::Matrix4d::Ones();
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num) {
        diff_se = Eigen::Matrix4d::Zero();
        for (int i = 0; i < N; i++) {
            diff_se += (Mean.inverse() * X[i]).log();
        }
        Mean *= ((1.0 / N)* diff_se).exp();
        count++;
    }

    // Covariance
    for (int i = 0; i < N; i++) {
        diff_se = (Mean.inverse() * X[i]).log();
        Eigen::VectorXd diff_vex(6);
        diff_vex << vex(diff_se.block<3, 3>(0, 0)), diff_se.block<3, 1>(0, 3);
        Cov += diff_vex * diff_vex.transpose();
    }
    Cov /= N;
}

Eigen::Matrix3d randomRotationMatrix(std::default_random_engine &generator) {
    std::normal_distribution<double> distribution(0.0, 1.0);
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m(i, j) = distribution(generator);
        }
    }

    Eigen::HouseholderQR<Eigen::Matrix3d> qr(m);
    Eigen::Matrix3d R = qr.householderQ();
    return R;
}

int main() {
    int num_matrices = 2;
    std::vector<Eigen::Matrix4d> X(num_matrices);

    X[0] << 0.9363, -0.2751,  0.2183,  1.2020,
            0.2896,  0.9566, -0.0392, -0.1022,
            -0.1985,  0.0978,  0.9750,  0.3426,
            0.0,     0.0,     0.0,     1.0;

    X[1] << 0.9938, -0.0975,  0.0599, -0.2246,
            0.0975,  0.9951, -0.0273,  0.1088,
            -0.0603,  0.0250,  0.9981,  0.4839,
            0.0,     0.0,     0.0,     1.0;

    // Call meanCov function
    Eigen::Matrix4d Mean;
    Eigen::Matrix<double, 6, 6> Cov;
    meanCov(X, Mean, Cov);

    // Display results
    std::cout << "Mean of transformation matrices:" << std::endl;
    std::cout << Mean << std::endl;
    std::cout << "Covariance of transformation matrices:" << std::endl;
    std::cout << Cov << std::endl;

    return 0;
}
