/*
DESCRIPTION:

The program defines a function skewLog that computes the matrix logarithm 
of a 3x3 rotation matrix R. The function first computes the angle of rotation 
theta using the trace of the matrix, then uses different cases depending on 
the value of theta to compute the corresponding skew-symmetric matrix w_hat. 
Finally, the computed w_hat matrix is returned. 
*/

#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>

Eigen::Matrix3d skewLog(Eigen::Matrix3d R) {
    Eigen::Matrix3d w_hat;
    double val = (R.trace() - 1.0) / 2.0;
    if (val > 1.0) {
        val = 1.0;
    } else if (val < -1.0) {
        val = -1.0;
    }
    double theta = std::acos(val);
    if (theta == 0) {
        w_hat.setZero();
    } else if (std::abs(M_PI - theta) < 1e-6) {
        Eigen::Matrix3d M = (R - Eigen::Matrix3d::Identity()) / 2.0;
        double m1 = M(0, 0);
        double m2 = M(1, 1);
        double m3 = M(2, 2);
        w_hat << 0, -std::sqrt((m3 - m1 - m2) / 2.0), std::sqrt((m2 - m1 - m3) / 2.0),
        std::sqrt((m3 - m1 - m2) / 2.0), 0, -std::sqrt((m1 - m2 - m3) / 2.0),
        -std::sqrt((m2 - m1 - m3) / 2.0), std::sqrt((m1 - m2 - m3) / 2.0), 0;
        w_hat *= theta;
    } else {
        w_hat = (R - R.transpose()) / (2.0 * std::sin(theta)) * theta;
    }
    return w_hat;
}