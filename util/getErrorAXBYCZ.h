/*
This defines a function getErrorAXBYCZ that takes in six Matrix4d objects 
representing rotation and translation transformations for three coordinate 
systems (X, Y, and Z). The function calculates the difference between the 
trace of each transformed matrix and the corresponding actual matrix and 
takes the arc-cosine of half of this difference. This is done for each of 
the three coordinate systems and the resulting values are stored in the 
first three elements of a Vector3d object. The function also calculates 
the Euclidean distance between the translation vectors of the transformed 
and actual matrices for each of the three coordinate systems, and these 
values are stored in the last three elements of the Vector3d object. 
The function returns this Vector3d object.

This code appears to be part of a larger program that likely involves 
transforming or registering data from multiple coordinate systems. 
The getErrorAXBYCZ function seems to calculate the differences between 
the transformed and actual coordinate systems and provides a way to 
quantify the error in the transformation. This information could be 
useful for evaluating the accuracy of the transformation and potentially 
refining the transformation to improve accuracy.
*/

#ifndef GETERRORAXBYCZ_H
#define GETERRORAXBYCZ_H

#include <eigen3/Eigen/Dense>
#include <cmath>

Eigen::Vector3d getErrorAXBYCZ(const Eigen::Matrix4d& X_f,
                               const Eigen::Matrix4d& Y_f,
                               const Eigen::Matrix4d& Z_f,
                               const Eigen::Matrix4d& XActual,
                               const Eigen::Matrix4d& YActual,
                               const Eigen::Matrix4d& ZActual) {

    Eigen::Vector3d xyzError;

    xyzError(0) = std::acos(0.5 * (X_f.block(0, 0, 3, 3).trace() - XActual.block(0, 0, 3, 3).trace()));
    xyzError(1) = std::acos(0.5 * (Y_f.block(0, 0, 3, 3).trace() - YActual.block(0, 0, 3, 3).trace()));
    xyzError(2) = std::acos(0.5 * (Z_f.block(0, 0, 3, 3).trace() - ZActual.block(0, 0, 3, 3).trace()));

    xyzError(3) = (X_f.block(0, 3, 3, 1) - XActual.block(0, 3, 3, 1)).norm();
    xyzError(4) = (Y_f.block(0, 3, 3, 1) - YActual.block(0, 3, 3, 1)).norm();
    xyzError(5) = (Z_f.block(0, 3, 3, 1) - ZActual.block(0, 3, 3, 1)).norm();

    return xyzError;
}

#endif