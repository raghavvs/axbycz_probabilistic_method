#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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
    Mean = expm(1.0 / N * sum_se);                  // use Eigen::MatrixBase::exp() for matrix exponential evaluation

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