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
#include <random>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

void meanCov(const Eigen::MatrixXd &X, int N, Eigen::MatrixXd &Mean, Eigen::MatrixXd &Cov) {
    std::vector<Eigen::MatrixXd> Y(N, X);
    Mean = Eigen::Matrix4d::Identity();
    Cov = Eigen::Matrix<double, 6, 6>::Zero();

    // Initial approximation of Mean
    Eigen::Matrix4d sum_se = Eigen::Matrix4d::Zero();
    for (int i = 0; i < N; i++) {
        sum_se += Y[i].log();
    }
    Mean = (1.0 / N * sum_se).exp();

    // Iterative process to calculate the true Mean
    Eigen::Matrix4d diff_se = Eigen::Matrix4d::Ones();
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num) {
        diff_se = Eigen::Matrix4d::Zero();
        for (int i = 0; i < N; i++) {
            diff_se += (Mean.inverse() * Y[i]).log();
        }
        Mean *= (1.0 / N * diff_se).exp();
        count++;
    }

    // Covariance
    for (int i = 0; i < N; i++) {
        diff_se = (Mean.inverse() * Y[i]).log();
        Eigen::VectorXd diff_vex(6);
        diff_vex << Eigen::Map<Eigen::Vector3d>(diff_se.block<3,3>(0,0).data()), diff_se.block<3,1>(0,3);
        Cov += diff_vex * diff_vex.transpose();
    }
    Cov /= N;
}

int main()
{
    int N = 1;
    Eigen::Matrix4d X = Eigen::Matrix4d::Random();
    std::cout << "X: " << X << std::endl;

    // Declare variables to store the output mean and covariance
    Eigen::MatrixXd Mean;
    Eigen::MatrixXd Cov;

    // Call the meanCov function with the input and output arguments
    meanCov(X, N,Mean,Cov);

    // Print the output mean and covariance
    std::cout << "The output mean is: " << std::endl;
    std::cout << Mean << std::endl;

    std::cout <<"The output covariance is:"<<std::endl;

    std::cout<<Cov<<std::endl;

    return 0;
}