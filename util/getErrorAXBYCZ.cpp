/*
DESCRIPTION:

The code defines a function getErrorAXBYCZ that takes in six Matrix4d objects 
(X_f, Y_f, Z_f, XActual, YActual, ZActual) and calculates the error between 
their orientations and translations using the cosine inverse and norm functions 
from the cmath library. The calculated error values are then stored in a 
Vector3d object (xyzError) and returned by the function.

Specifically, for each input matrix, the function calculates the orientation 
error using the difference in trace values between the matrix and the corresponding 
actual matrix. The translation error is calculated using the norm of the difference 
between the translation vector of the matrix and the corresponding actual matrix. 
The calculated error values are stored in the xyzError vector and returned.

std::acos to calculate the rotational errors
*/

#include <iostream>
#include <Eigen/Dense>
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

// TEST

int main() {
    Eigen::Matrix4d X_f, Y_f, Z_f, XActual, YActual, ZActual;
    X_f << 1, 0, 0, 2,
           0, 1, 0, 3,
           0, 0, 1, 4,
           0, 0, 0, 1;
    Y_f << 1, 0, 0, 2,
           0, 1, 0, 3,
           0, 0, 1, 4,
           0, 0, 0, 1;
    Z_f << 1, 0, 0, 2,
           0, 1, 0, 3,
           0, 0, 1, 4,
           0, 0, 0, 1;
    XActual << 1, 0, 0, 2,
               0, 1, 0, 3,
               0, 0, 1, 4,
               0, 0, 0, 1;
    YActual << 1, 0, 0, 2,
               0, 1, 0, 3,
               0, 0, 1, 4,
               0, 0, 0, 1;
    ZActual << 1, 0, 0, 2,
               0, 1, 0, 3,
               0, 0, 1, 4,
               0, 0, 0, 1;

    Eigen::Vector3d expected = Eigen::Vector3d::Zero();
    Eigen::Vector3d actual = getErrorAXBYCZ(X_f, Y_f, Z_f, XActual, YActual, ZActual);

    if (expected.isApprox(actual, 1e-5)) {
        std::cout << "Test passed" << std::endl;
    } else {
        std::cout << "Test failed" << std::endl;
    }

    return 0;
}