#include <eigen3/Eigen/Dense>
#include <cmath>

Eigen::Vector3d getErrorAXBYCZ(const Eigen::Matrix4d& X_f, const Eigen::Matrix4d& Y_f, const Eigen::Matrix4d& Z_f,
                               const Eigen::Matrix4d& XActual, const Eigen::Matrix4d& YActual, const Eigen::Matrix4d& ZActual) {

    Eigen::Vector3d xyzError;

    xyzError(0) = std::acos(0.5 * (X_f.block(0, 0, 3, 3).trace() - XActual.block(0, 0, 3, 3).trace()));
    xyzError(1) = std::acos(0.5 * (Y_f.block(0, 0, 3, 3).trace() - YActual.block(0, 0, 3, 3).trace()));
    xyzError(2) = std::acos(0.5 * (Z_f.block(0, 0, 3, 3).trace() - ZActual.block(0, 0, 3, 3).trace()));

    xyzError(3) = (X_f.block(0, 3, 3, 1) - XActual.block(0, 3, 3, 1)).norm();
    xyzError(4) = (Y_f.block(0, 3, 3, 1) - YActual.block(0, 3, 3, 1)).norm();
    xyzError(5) = (Z_f.block(0, 3, 3, 1) - ZActual.block(0, 3, 3, 1)).norm();

    return xyzError;
}

/*
std::acos to calculate the rotational errors, and norm() function to calculate the translational errors.
*/