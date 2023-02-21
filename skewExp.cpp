#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

Matrix3d skew(Vector3d v)
{
    Matrix3d m;
    m << 0, -v(2), v(1),
         v(2), 0, -v(0),
         -v(1), v(0), 0;
    return m;
}

Matrix3d skewExp(Vector3d s, double theta = 1)
{
    Matrix3d g;
    g.setIdentity();
    Matrix3d ss = skew(s);
    g += ss * sin(theta) + ss * ss * (1 - cos(theta));
    return g;
}