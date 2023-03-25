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
#include <ctime>
#include <random>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

void meanCov(const std::vector<Eigen::Matrix4d> &X, int N, std::vector<Eigen::MatrixXd> &Mean,
             std::vector<Eigen::MatrixXd> &Cov) {
    Mean[0] = Eigen::MatrixXd::Identity(4, 4);
    Cov[0] = Eigen::MatrixXd::Zero(6, 6);

    // Initial approximation of Mean
    Eigen::MatrixXd sum_se = Eigen::MatrixXd::Zero(4, 4);
    for (int i = 0; i < N; i++) {
        sum_se += X[i].log();
    }
    Mean[0] = (1.0 / N * sum_se).exp();

    // Iterative process to calculate the true Mean
    Eigen::MatrixXd diff_se = Eigen::MatrixXd::Ones(4, 4);
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num) {
        diff_se = Eigen::MatrixXd::Zero(4, 4);
        for (int i = 0; i < N; i++) {
            diff_se += (Mean[0].inverse() * X[i]).log();
        }
        Mean[0] *= (1.0 / N * diff_se).exp();
        count++;
    }
    // Covariance
    for (int i = 0; i < N; i++) {
        diff_se = (Mean[0].inverse() * X[i]).log();
        Eigen::VectorXd diff_vex(6);
        diff_vex << Eigen::Map<Eigen::Vector3d>(diff_se.block<3,3>(0,0).data()), diff_se.block<3,1>(0,3);
        Cov[0] += diff_vex * diff_vex.transpose();
    }
    Cov[0] /= N;
}

int main() {
    // Example data
    std::vector<Eigen::Matrix4d> X(3);
    X[0] << 1.1, 2.2, 3.3, 4.4,
            5.5, 6.6, 7.7, 8.8,
            9.9, 10.10, 11.11, 12.12,
            13.13, 14.14, 15.15, 16.16;
    X[1] << -1.1, -2.2, -3.3, -4.4,
            -5.5, -6.6, -7.7, -8.8,
            -9.9, -10.10, -11.11, -12.12,
            -13.13, -14.14,-15.15,-16.16;
    X[2] << 2,-2,-2,-2,
            -2,-2,-2,-2,
            -2,-2,-2,-2,
            -2,-2,-2,-2;

    int N = X.size();

    std::vector<Eigen::MatrixXd> Mean(N);
    std::vector<Eigen::MatrixXd> Cov(N);

    meanCov(X,N ,Mean,Cov);

    std::cout << "Mean:\n" << Mean[2] << "\n\n";
    std::cout << "Covariance:\n" << Cov[2] << "\n\n";

    return 0;
}