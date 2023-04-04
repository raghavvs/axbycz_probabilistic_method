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

#include <iostream>
#include <eigen3/Eigen/Dense>
#include "rotError.h"
#include "tranError.h"

Eigen::VectorXd getErrorAXBYCZ(const Eigen::Matrix4d& X_f, 
                               const Eigen::Matrix4d& Y_f, 
                               const Eigen::Matrix4d& Z_f,
                               const Eigen::Matrix4d& XActual, 
                               const Eigen::Matrix4d& YActual, 
                               const Eigen::Matrix4d& ZActual) {
    Eigen::VectorXd xyzError(6);
    
    xyzError(0) = rotError(X_f, XActual);
    xyzError(1) = rotError(Y_f, YActual);
    xyzError(2) = rotError(Z_f, ZActual);

    xyzError(3) = tranError(X_f, XActual);
    xyzError(4) = tranError(Y_f, YActual);
    xyzError(5) = tranError(Z_f, ZActual);

    return xyzError;
}

// std::acos to calculate the rotational errors, and norm() function to calculate the translational errors.

int main() {
    Eigen::Matrix4d X_f = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d Y_f = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d Z_f = Eigen::Matrix4d::Identity();

    Eigen::Matrix4d XActual = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d YActual = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d ZActual = Eigen::Matrix4d::Identity();

    // Introduce some differences
    XActual(0, 0) += 0.1;
    XActual(0, 1) += 0.8;
    YActual(1, 1) += 0.2;
    YActual(2, 1) += 0.7;
    ZActual(2, 2) += 0.3;

    XActual(0, 3) += 0.4;
    YActual(1, 3) += 0.5;
    ZActual(2, 3) += 0.6;

    Eigen::VectorXd xyzError = getErrorAXBYCZ(X_f, Y_f, Z_f, XActual, YActual, ZActual);

    std::cout << "Error vector: \n" << xyzError << std::endl;

    return 0;
}