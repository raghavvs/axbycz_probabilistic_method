#include <eigen3/Eigen/Dense>

double rotError(Eigen::Matrix4d X1, Eigen::Matrix4d X2)
{
    Eigen::Matrix3d R1 = X1.block<3, 3>(0, 0);
    Eigen::Matrix3d R2 = X2.block<3, 3>(0, 0);
    Eigen::Matrix3d R12 = R1.transpose() * R2;
    Eigen::AngleAxisd aa(R12);
    return aa.angle();
}
