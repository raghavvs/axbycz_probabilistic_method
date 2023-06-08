//
// Created by Raghavendra N S on 6/7/23.
//

#ifndef PARAMEXTRACT_H
#define PARAMEXTRACT_H

#include <Eigen/Dense>
#include <iostream>

Eigen::MatrixXd so3Vec(const Eigen::MatrixXd& X)
{
    if (X.cols() == 3) // If input is skew-sym change to vector
    {
        Eigen::VectorXd g(3);
        g << -X(1,2), X(0,2), -X(0,1);
        return g;
    }
    else // If input is vector change to skew-sym
    {
        Eigen::Matrix3d g;
        g << 0, -X(2), X(1),
                X(2), 0, -X(0),
                -X(1), X(0), 0;
        return g;
    }
}

void paramExtract(double& theta,
                  Eigen::MatrixXd& N,
                  double& d,
                  Eigen::VectorXd& p,
                  const Eigen::MatrixXd& X) {
    // Extract theta
    theta = std::acos((X.block<3,3>(0,0).trace() - 1) / 2);

    // Extract N
    N = (X.block<3,3>(0,0) - X.block<3,3>(0,0).transpose()) / (2 * std::sin(theta));

    // Extract d
    Eigen::Vector3d n = so3Vec(N);
    d = X.block<3,1>(0,3).dot(n);

    // Extract p
    Eigen::Vector3d u(3);
    u << -n(1), n(0), 0;
    u = (1 / std::sqrt(n(0)*n(0) + n(1)*n(1))) * u;
    Eigen::MatrixXd A(2,2);
    A << 1 - std::cos(theta), std::sin(theta),
            -std::sin(theta), 1 - std::cos(theta);
    Eigen::VectorXd b(2);
    Eigen::Vector3d u3 = u;
    Eigen::Vector3d crossProduct = n.cross(u);
    b << X.block<3,1>(0,3).dot(u3), X.block<3,1>(0,3).dot(crossProduct);
    Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);

    p = c(0) * u + c(1) * n.cross(u);
}

#endif
