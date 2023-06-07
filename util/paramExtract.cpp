/*
DESCRIPTION:

*/

#include "paramExtract.h"

void paramExtract(double& theta,
                  const& Eigen::MatrixXd N,
                  double& d
                  Eigen::VectorXd& p,
                  const& Eigen::MatrixXd X) {
    // Extract theta
    theta = std::acos((X.block<3,3>(0,0).trace() - 1) / 2);

    // Extract N
    N = (X.block<3,3>(0,0) - X.block<3,3>(0,0).transpose()) / (2 * std::sin(theta));

    // Extract d
    d = X.block<3,1>(0,4).dot(so3_vec(N));

    // Extract p
    Eigen::VectorXd n = so3Vec(N);
    Eigen::VectorXd u = (1 / std::sqrt(n(0)*n(0) + n(1)*n(1))) * Eigen::VectorXd{-n(1), n(0), 0};
    Eigen::MatrixXd A(2,2);
    A << 1 - std::cos(theta), std::sin(theta),
            -std::sin(theta), 1 - std::cos(theta);
    Eigen::VectorXd b(2);
    b << X.block<3,1>(0,4).dot(u), X.block<3,1>(0,4).dot(n.cross(u));
    Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);

    p = c(0) * u + c(1) * n.cross(u);
}