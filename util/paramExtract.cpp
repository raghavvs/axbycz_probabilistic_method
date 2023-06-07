/*
DESCRIPTION:


*/

#include <Eigen/Dense>
#include "so3Vec.h"

// Assuming the custom functions are defined elsewhere
Eigen::VectorXd so3Cec(const Eigen::MatrixXd& N);

void param_extract(const Eigen::MatrixXd& X) {
    // Extract theta
    double theta = std::acos((X.block<3,3>(0,0).trace() - 1) / 2);

    // Extract N
    Eigen::MatrixXd N = (X.block<3,3>(0,0) - X.block<3,3>(0,0).transpose()) / (2 * std::sin(theta));

    // Extract d
    double d = X.block<3,1>(0,4).dot(so3_vec(N));

    // Extract p
    Eigen::VectorXd n = so3_vec(N);
    Eigen::VectorXd u = (1 / std::sqrt(n(0)*n(0) + n(1)*n(1))) * Eigen::VectorXd{-n(1), n(0), 0};
    Eigen::MatrixXd A(2,2);
    A << 1 - std::cos(theta), std::sin(theta),
            -std::sin(theta), 1 - std::cos(theta);
    Eigen::VectorXd b(2);
    b << X.block<3,1>(0,4).dot(u), X.block<3,1>(0,4).dot(n.cross(u));
    Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);  // Solve the linear system

    Eigen::VectorXd p = c(0) * u + c(1) * n.cross(u);

    // The variables theta, N, d, and p now hold the results
}
