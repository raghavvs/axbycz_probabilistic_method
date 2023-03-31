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

void meanCov(const std::vector<Eigen::Matrix4d> &X,
             int N,
             std::vector<Eigen::MatrixXd> &Mean,
             std::vector<Eigen::MatrixXd> &Cov) {
    for(int i = 0; i < N; i++){
        Mean.push_back(Eigen::Matrix4d::Identity());
        Cov.push_back(Eigen::Matrix<double, 6, 6>::Zero());
    }

    // Initial approximation of Mean
    Eigen::Matrix4d sum_se = Eigen::Matrix4d::Zero();
    for (int i = 0; i < N; i++) {
        sum_se += X[i].log();
        Mean[i] = ((1.0 / N) * sum_se).exp();
    }

    // Iterative process to calculate the true Mean
    Eigen::Matrix4d diff_se = Eigen::Matrix4d::Ones();
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num) {
        diff_se = Eigen::Matrix4d::Zero();
        for (int i = 0; i < N; i++) {
            diff_se += (Mean[i].inverse() * X[i]).log();
            Mean[i] *= ((1.0 / N)* diff_se).exp();
        }
        count++;
    }

    // Covariance
    for (int i = 0; i < N; i++) {
        diff_se = (Mean[i].inverse() * X[i]).log();
        Eigen::VectorXd diff_vex(6);
        diff_vex << Eigen::Map<Eigen::Vector3d>(diff_se.block<3,3>(0,0).data()), diff_se.block<3,1>(0,3);
        Cov[i] += diff_vex * diff_vex.transpose();
        Cov[i] /= N;
    }
}

int main() {
    std::vector<Eigen::Matrix4d> A1;
    std::vector<Eigen::Matrix4d> B1;
    std::vector<Eigen::Matrix4d> C1;
    std::vector<Eigen::Matrix4d> A2;
    std::vector<Eigen::Matrix4d> B2;
    std::vector<Eigen::Matrix4d> C2;

    for (int i = 0; i < 10; ++i) {
        A1.emplace_back(Eigen::Matrix4d::Random());
        B1.emplace_back(Eigen::Matrix4d::Random());
        C1.emplace_back(Eigen::Matrix4d::Random());
        A2.emplace_back(Eigen::Matrix4d::Random());
        B2.emplace_back(Eigen::Matrix4d::Random());
        C2.emplace_back(Eigen::Matrix4d::Random());
    }

    int Ni = A1.size();
    int Nj = C2.size();

    std::vector<Eigen::MatrixXd> A1_m, B1_m, C1_m, SigA1, SigB1, SigC1;
    meanCov(A1, Ni, A1_m, SigA1);
    meanCov(B1, Ni, B1_m, SigB1);
    meanCov(C1, Ni, C1_m, SigC1);

    std::cout << "Mean:\n" << A1_m[2] << std::endl;
    std::cout << "Covariance:\n" << SigA1[2] << std::endl;

    return 0;
}