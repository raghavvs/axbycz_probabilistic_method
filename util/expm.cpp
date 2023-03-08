/*
DESCRIPTION:

A function to evaluate matrix exponential.
It uses the scaling and squaring method.
Process:
    Compute the scaling factor s such that ||A/s|| <= 1/2, where ||.|| is a matrix norm (e.g., the Frobenius norm).
    Compute the matrix B = A/s.
    Compute the matrix exponential of B using the power series expansion.
    Compute the matrix exponential of A as exp(A) = (exp(s) * exp(B))^2.
n term refers to number of terms in the power expansion.
*/

#include <cmath>
#include <eigen3/Eigen/Dense>

Eigen::MatrixXd expm(const Eigen::MatrixXd& A, int nterms = 20) {
    double s = 0.5;
    double normA = A.norm();
    while (normA / s > 1.0) {
        s *= 2.0;
    }
    Eigen::MatrixXd B = A / s;
    Eigen::MatrixXd X = Eigen::MatrixXd::Identity();
    Eigen::MatrixXd C = Eigen::MatrixXd::Identity();
    for (int k = 1; k <= nterms; k++) {
        X *= B / k;
        C += X;
    }
    for (int i = 0; i < static_cast<int>(s); i++) {
        C = C * C;
    }
    return C;
}