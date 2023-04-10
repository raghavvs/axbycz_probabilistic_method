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
*/

#ifndef MEANCOV_H
#define MEANCOV_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <unsupported/Eigen/MatrixFunctions>

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
        Eigen::Vector3d diff_se_block_vec = Eigen::Vector3d::Map(diff_se.block<3, 3>(0, 0).data());
        diff_vex << diff_se_block_vec, diff_se.block<3, 1>(0, 3);
        Cov += diff_vex * diff_vex.transpose();
    }
    Cov /= N;
}

#endif