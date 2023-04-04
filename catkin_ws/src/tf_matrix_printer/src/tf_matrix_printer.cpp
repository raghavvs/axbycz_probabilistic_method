//
// Created by raghav on 4/3/23.
//

#include <iostream>
#include <tf/transform_datatypes.h>

int main() {
    // Assume that tf_echo returns the following translation vector and quaternion
    double x = 1.0;
    double y = 2.0;
    double z = 3.0;
    double qx = 0.0;
    double qy = 0.0;
    double qz = 0.0;
    double qw = 1.0;

    // Convert the translation vector into a tf::Vector3 object
    tf::Vector3 translation(x, y, z);

    // Convert the quaternion into a tf::Quaternion object
    tf::Quaternion rotation(qx, qy, qz, qw);

    // Create a tf::Transform object representing the transformation between the two frames
    tf::Transform transform;
    transform.setOrigin(translation);
    transform.setRotation(rotation);

    // Convert the rotation part of the transform to a 3x3 matrix
    tf::Matrix3x3 rotation_matrix(transform.getRotation());

    // Print the transformation matrix
    std::cout << "Transformation matrix: " << std::endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << rotation_matrix[i][j] << " ";
        }
        // Print the translation part
        std::cout << translation[i] << std::endl;
    }
    // Print the last row (for a homogeneous transformation matrix)
    std::cout << "0 0 0 1" << std::endl;

    return 0;
}