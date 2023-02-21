#include <iostream>
#include <eigen3/Eigen/Dense>

double tranError(Eigen::MatrixXd X1, Eigen::MatrixXd X2) {
    Eigen::Vector3d p1 = X1.block<3,1>(0,3);
    Eigen::Vector3d p2 = X2.block<3,1>(0,3);

    return (p1 - p2).norm();
}