/*
DESCRIPTION:

This program defines several functions for calculating the mean 
and covariance of a set of 4x4 matrices. The logm function takes 
a 4x4 matrix and calculates its matrix logarithm, while the vex 
function takes a 3x3 matrix and returns its vector of exteriorization. 
The meanCov function takes an array of 3x3 matrices and its size 
N and calculates the mean and covariance of the logarithms of those 
matrices. It does this by first taking the average of the logarithms 
using the expm function, then iteratively refining this average until 
convergence using the logm and vex functions. Finally, it calculates 
the covariance by taking the vector of differences between each logarithm 
and the mean logarithm and computing their outer product.
*/

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <expm.h>
#include <iostream>

Eigen::Matrix4d logm(Eigen::Matrix4d A) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> es(A);
    Eigen::Matrix4d D = es.eigenvalues().asDiagonal();
    Eigen::Matrix4d V = es.eigenvectors();
    Eigen::Matrix4d S = V * (D.array().log().matrix().asDiagonal()) * V.inverse();
    return S;
}

Eigen::VectorXd vex(Eigen::Matrix3d S) {
    Eigen::VectorXd w(3);
    w(0) = S(2, 1);
    w(1) = S(0, 2);
    w(2) = S(1, 0);
    return w;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> meanCov(const Eigen::MatrixXd& X, Eigen::VectorXd& Mean, Eigen::MatrixXd& Cov)
{
    int N = X.size() / 16; // size of third dimension of X
    Mean = Eigen::Matrix4d::Identity();
    Cov = Eigen::MatrixXd::Zero(6, 6);

    // Initial approximation of Mean
    Eigen::Matrix4d sum_se = Eigen::Matrix4d::Zero();
    for (int i = 0; i < N; i++)
    {
        Eigen::Matrix4d Xi = Eigen::Map<const Eigen::Matrix4d>(X.data() + i * 16, 4, 4);
        sum_se += logm(Xi);
    }
    Mean = expm(1.0 / N * sum_se);

    // Iterative process to calculate the true Mean
    Eigen::Matrix4d diff_se = Eigen::Matrix4d::Ones();
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num)
    {
        diff_se = Eigen::Matrix4d::Zero();
        for (int i = 0; i < N; i++)
        {
            Eigen::Matrix4d Xi = Eigen::Map<const Eigen::Matrix4d>(X.data() + i * 16, 4, 4);
            diff_se += logm(Mean.inverse() * Xi);
        }
        Mean *= expm(1.0 / N * diff_se);
        count++;
    }

    Eigen::VectorXd diff_vex;

    // Covariance
    for (int i = 0; i < N; i++)
    {
        Eigen::Matrix4d Xi = Eigen::Map<const Eigen::Matrix4d>(X.data() + i * 16, 4, 4);
        Eigen::Matrix4d diff_se = logm(Mean.inverse() * Xi);
        diff_vex << vex(diff_se.block<3, 3>(0, 0)), diff_se.block<3, 1>(0, 3);
        Cov += diff_vex * diff_vex.transpose();
    }
    Cov /= N;

    return std::make_pair(vex(logm(Mean)), Cov);
}

/*
vex() is a function that takes a skew-symmetric matrix 
and returns its corresponding 3D vector.
*/