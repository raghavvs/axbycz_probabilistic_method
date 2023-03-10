/*
DESCRIPTION:

The program defines two functions: skew() and skewExp(). The skew() 
function takes in a 3D vector and returns a 3x3 skew-symmetric matrix. 
The skewExp() function takes in a 3D vector s and an angle theta 
(default value is 1), and returns a 3x3 matrix calculated using the 
exponential map of the 3D vector. The matrix is constructed using the 
skew-symmetric matrix of s, and the rotation matrix calculated using 
the Rodrigues formula.
*/

#include <iostream>
#include <eigen3/Eigen/Dense>

Eigen::Matrix3d skew(Eigen::Vector3d v)
{
    Eigen::Matrix3d m;
    m << 0, -v(2), v(1),
         v(2), 0, -v(0),
         -v(1), v(0), 0;
    return m;
}

Eigen::Matrix3d skewExp(Eigen::Vector3d s, double theta = 1)
{
    Eigen::Matrix3d g;
    g.setIdentity();
    Eigen::Matrix3d ss = skew(s);
    g += ss * sin(theta) + ss * ss * (1 - cos(theta));
    return g;
}