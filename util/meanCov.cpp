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
//#include <expm.h>
#include <cmath>
#include <iostream>

Eigen::MatrixXd expm(const Eigen::MatrixXd& A, int nterms = 20) {
    double s = 0.5;
    double normA = A.norm();
    while (normA / s > 1.0) {
        s *= 2.0;
    }
    Eigen::MatrixXd B = A / s;
    Eigen::MatrixXd X = Eigen::MatrixXd::Identity(A.rows(), A.cols());
    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(A.rows(), A.cols());
    for (int k = 1; k <= nterms; k++) {
        X *= B / k;
        C += X;
    }
    for (int i = 0; i < static_cast<int>(s); i++) {
        C = C * C;
    }
    return C;
}

Eigen::MatrixXd logm(Eigen::MatrixXd A) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    Eigen::MatrixXd D = es.eigenvalues().asDiagonal();
    Eigen::MatrixXd V = es.eigenvectors();
    Eigen::MatrixXd S = V * (D.array().log().matrix().asDiagonal()) * V.inverse();
    return S;
}

Eigen::Vector3d vex(const Eigen::Matrix3d& S){
    Eigen::Vector3d w;
    w << S(2, 1), S(0, 2), S(1, 0);
    return w;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> meanCov(const Eigen::MatrixXd& X, Eigen::VectorXd& Mean, Eigen::MatrixXd& Cov)
{
    int N = X.rows() / 4; // size of third dimension of X
    Mean.setIdentity(4, 4);
    Cov.setZero(6, 6);

    // Initial approximation of Mean
    Eigen::MatrixXd sum_se = Eigen::MatrixXd::Zero(4, 4);
    for (int i = 0; i < N; i++)
    {
        Eigen::MatrixXd Xi = Eigen::Map<const Eigen::MatrixXd>(X.data() + i * 16, 4, 4);
        sum_se += logm(Xi);
    }
    Mean = expm(1.0 / N * sum_se);

    // Iterative process to calculate true Mean
    Eigen::MatrixXd diff_se = Eigen::MatrixXd::Ones(4, 4);
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num)
    {
        diff_se = Eigen::MatrixXd::Zero(4, 4);
        for (int i = 0; i < N; i++)
        {
            Eigen::MatrixXd Xi = Eigen::Map<const Eigen::MatrixXd>(X.data() + i * 16, 4, 4);
            diff_se += logm(Mean.inverse() * Xi);
        }
        Mean *= expm(1.0 / N * diff_se);
        count++;
    }

    Eigen::VectorXd diff_vex;

    // Covariance
    for (int i = 0; i < N; i++)
    {
        Eigen::MatrixXd Xi = Eigen::Map<const Eigen::MatrixXd>(X.data() + i * 16, 4, 4);
        Eigen::MatrixXd diff_se = logm(Mean.inverse() * Xi);
        diff_vex << vex(diff_se.block<3, 3>(0, 0).template cast<double>()), diff_se.block<3, 1>(0, 3);
        Cov = Cov + diff_vex * diff_vex.transpose();
    }
    Cov = Cov / N;

    return std::make_pair(vex(logm(Mean)), Cov);
}

/*
vex() is a function that takes a skew-symmetric matrix 
and returns its corresponding 3D vector.
*/

//TEST Case

/* int main()
{
    // Define some 4x4 matrices
    Eigen::MatrixXd X1;
    X1 << 1, 0, 0, 1,
          0, 1, 0, 2,
          0, 0, 1, 3,
          0, 0, 0, 1;
          
    Eigen::MatrixXd X2;
    X2 << -1, 0, 0, 4,
          0, -1, 0, 5,
          0, 0, -1, 6,
          0, 0, 0, 1;
          
    Eigen::MatrixXd X3;
    X3 << 0, 1, 0, 7,
          -1, 0, 0, 8,
          0, 0, 1, 9,
          0, 0, 0, 1;

    // Put the matrices in a vector
    Eigen::MatrixXd X(16, 3);
    Eigen::Map<Eigen::MatrixXd>(X.data(), 4, 4) = X1;
    Eigen::Map<Eigen::MatrixXd>(X.data() + 16, 4, 4) = X2;
    Eigen::Map<Eigen::MatrixXd>(X.data() + 32, 4, 4) = X3;

    // Compute the mean and covariance
    Eigen::VectorXd Mean(4);
    Eigen::MatrixXd Cov(6, 6);
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> result = meanCov(X, Mean, Cov);

    // Print the results
    std::cout << "Mean:\n" << result.first << "\n";
    std::cout << "Covariance:\n" << result.second << "\n";

    return 0;
} */