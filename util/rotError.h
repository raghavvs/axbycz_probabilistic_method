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

double rotError(Eigen::Matrix4d X1, Eigen::Matrix4d X2)
{
    Eigen::Matrix3d R1 = X1.block<3, 3>(0, 0);
    Eigen::Matrix3d R2 = X2.block<3, 3>(0, 0);
    Eigen::Matrix3d R12 = R1.transpose() * R2;
    Eigen::AngleAxisd aa(R12);
    return aa.angle();
}

// TEST CASE

/*
int main()
{
    Eigen::Matrix4d X1, X2;

    X1 << 0.7204, -0.0024, -0.6937, -0.7254,
          0.6936, -0.1371,  0.7070,  0.0161,
          -0.0092, -0.9905, -0.1375, -0.2605,
           0.0000,  0.0000,  0.0000,  1.0000;
    X2 << 0.7185, -0.0045, -0.6955, -0.7225,
          0.6953, -0.1441,  0.7039,  0.0098,
          -0.0065, -0.9895, -0.1448, -0.2614,
           0.0000,  0.0000,  0.0000,  1.0000;
    
    std::cout << "Rotation error: " << rotError(X1, X2) << std::endl;

    return 0;
}
*/