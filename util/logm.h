#ifndef LOGM_H
#define LOGM_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <cmath>
#include <iostream>

Eigen::MatrixXd logm(const Eigen::MatrixXd& A)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();

    Eigen::MatrixXd D = svd.singularValues();

    for (int i = 0; i < D.rows(); ++i)
    {
        if (D(i) < 0)
        {
            D(i) = 0;
        }
        else
        {
            D(i) = std::log(D(i));
        }
    }

    return U * D.asDiagonal() * V.transpose();
}

#endif