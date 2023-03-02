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

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <expm.h>

using namespace Eigen;
using namespace std;

Matrix4d logm(Matrix4d A) {
    SelfAdjointEigenSolver<Matrix4d> es(A);
    Matrix4d D = es.eigenvalues().asDiagonal();
    Matrix4d V = es.eigenvectors();
    Matrix4d S = V * (D.array().log().matrix().asDiagonal()) * V.inverse();
    return S;
}

VectorXd vex(Matrix3d S) {
    VectorXd w(3);
    w(0) = S(2, 1);
    w(1) = S(0, 2);
    w(2) = S(1, 0);
    return w;
}

void meanCov(Matrix3d* X, int N, Matrix4d& Mean, MatrixXd& Cov) {
    Mean = Matrix4d::Identity();

    Matrix4d sum_se = Matrix4d::Zero();
    for (int i = 0; i < N; i++) {
        sum_se += logm(X[i]);
    }
    Mean = expm(1.0 / N * sum_se);                  

    Matrix4d diff_se = Matrix4d::Ones();
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num) {
        diff_se = Matrix4d::Zero();
        for (int i = 0; i < N; i++) {
            diff_se += logm(Mean.inverse() * X[i]);
        }
        Mean = Mean * expm(1.0 / N * diff_se);
        count++;
    }

    Cov = MatrixXd::Zero(6, 6);
    for (int i = 0; i < N; i++) {
        Matrix4d diff_se = logm(Mean.inverse() * X[i]);
        VectorXd diff_vex = vex(diff_se.block<3, 3>(0, 0));
        diff_vex.conservativeResize(6);
        diff_vex.tail(3) = diff_se.block<3, 1>(0, 3);
        Cov += diff_vex * diff_vex.transpose();
    }
    Cov /= N;
}

/*
vex() is a function that takes a skew-symmetric matrix and returns its corresponding 3D vector.
*/