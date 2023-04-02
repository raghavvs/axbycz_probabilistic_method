/*
Description:

The code above defines a function named rotError that computes 
the rotation error between two 4x4 transformation matrices X1 
and X2. The rotation error is defined as the angle needed to 
rotate R1, the rotation matrix in X1, to align with R2, the 
rotation matrix in X2. To compute the rotation error, the 
function first extracts the rotation matrices R1 and R2 from 
the transformation matrices X1 and X2, respectively. It then 
computes the rotation matrix R12 that maps R1 to R2. Finally, 
it computes the angle-axis representation of R12 using the 
Eigen::AngleAxisd class and returns the angle component of 
the resulting object.

The Eigen::AngleAxisd class represents a rotation as an angle 
of rotation about a given axis in 3D space. Therefore, the 
rotError function returns the angle component of the rotation 
between the two matrices in radians.
*/

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>

Eigen::Matrix3d skewlog(const Eigen::Matrix3d& R) {
    double theta = acos((R.trace() - 1) / 2);
    Eigen::Matrix3d omega = (R - R.transpose()) / (2 * sin(theta));

    return theta * omega;
}

Eigen::Vector3d so3_vec(const Eigen::Matrix3d& omega) {
    return Eigen::Vector3d(omega(2, 1), omega(0, 2), omega(1, 0));
}

double rotError(const Eigen::Matrix4d& X1, const Eigen::Matrix4d& X2) {
    Eigen::Matrix3d R1 = X1.block<3, 3>(0, 0);
    Eigen::Matrix3d R2 = X2.block<3, 3>(0, 0);
    Eigen::Matrix3d R12 = R1.transpose() * R2;

    Eigen::Vector3d err_vec = so3_vec(skewlog(R12));

    return err_vec.norm();
}

// TEST CASE

int main() {
    // Create two example transformation matrices
    Eigen::Matrix4d X1, X2;

    X1 << 1, 0, 0, 1,
          0, 1, 0, 2,
          0, 0, 1, 3,
          0, 0, 0, 1;

    X2 << 0.866, -0.5, 0, 1.5,
          0.5, 0.866, 0, 2.5,
          0, 0, 1, 3.5,
          0, 0, 0, 1;

    // Call the rotError function to compute the rotational error
    double error = rotError(X1, X2);

    // Print the result
    std::cout << "Rotational error: " << error << std::endl;

    return 0;
}
