/*
DESCRIPTION:

The provided code implements a set of functions that solve a robotics problem 
involving the transformation matrices of multiple coordinate frames. Specifically, 
the functions solve for the transformations between three coordinate frames 
(A, B, and C) given the transformations between A and B and between A and C. 
The functions use the Eigen library to perform matrix operations such as inversion 
and SVD decomposition. The main function (axbyczProb1) calls the other two functions 
(batchSolveXY and randSE3) to generate a set of random transformations and iteratively 
select those that satisfy certain constraints, in order to estimate the desired transformations.
*/

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <rotError.h>

Eigen::Matrix4d invSE3(const Eigen::Matrix4d& g)
{
    Eigen::Matrix3d R = g.block<3, 3>(0, 0);
    Eigen::Vector3d t = g.block<3, 1>(0, 3);
    Eigen::Matrix4d gInv = Eigen::Matrix4d::Identity();
    gInv.block<3, 3>(0, 0) = R.transpose();
    gInv.block<3, 1>(0, 3) = -R.transpose() * t;
    return gInv;
}

Eigen::Matrix4d randSE3()
{
    Eigen::VectorXd u = Eigen::VectorXd::Random(6);
    Eigen::Matrix4d g = Eigen::Matrix4d::Identity();
    g.topLeftCorner(3, 3) = Eigen::AngleAxisd(u.tail<3>().norm(), u.tail<3>().normalized()).toRotationMatrix();
    g.block<3, 1>(0, 3) = 5 * u.head<3>();
    return g;
}

void batchSolveXY(Eigen::Matrix4d A, Eigen::Matrix4d B, Eigen::Matrix4d& X, Eigen::Matrix4d& Y, Eigen::MatrixXd& MeanA, Eigen::MatrixXd& MeanB, const double nstd1 = 2, const double nstd2 = 2)
{
    Eigen::MatrixXd g(4, 4);
    Eigen::MatrixXd samples(4, 4);
    samples.setZero();
    MeanA.setZero();
    MeanB.setZero();
    X.setZero();
    Y.setZero();

    int num_samples = 1000;
    int good_samples = 0;

    // sample rotations
    while (good_samples < num_samples)
    {
        g = randSE3();

        // check if g is within nstd1 of A
        if (rotError(g, A) < nstd1)
        {
            // sample translations
            for (int i = 0; i < 5; i++)
            {
                g.block<3, 1>(0, 3) = 5 * Eigen::Vector3d::Random();

                // check if g is within nstd2 of B
                if (rotError(g, B) < nstd2)
                {
                    samples += g;
                    MeanA += invSE3(g) * A;
                    MeanB += invSE3(g) * B;
                    good_samples++;
                    break;
                }
            }
        }
    }

    MeanA /= num_samples;
    MeanB /= num_samples;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(MeanB, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd Sigma(4, 4);
    Sigma.setZero();
    Sigma.topLeftCorner(3, 3) = svd.singularValues().asDiagonal();
    Sigma(3, 3) = 1;
    X = MeanA * svd.matrixV() * Sigma * svd.matrixU().transpose();
    Y = svd.matrixU() * Sigma * svd.matrixV().transpose();
}

void axbyczProb1(Eigen::Matrix4d A1, Eigen::Matrix4d B1, Eigen::Matrix4d C1, 
                Eigen::Matrix4d A2, Eigen::Matrix4d B2, Eigen::Matrix4d C2, 
                Eigen::Matrix4d& X_final, Eigen::Matrix4d& Y_final, Eigen::Matrix4d& Z_final, 
                const double opt, const double nstd1, const double nstd2)
{
    Eigen::Matrix4d X_temp, Y_temp, Z_temp;
    Eigen::Vector4d b1, b2;
    b1 << B1.col(3);
    b2 << B2.col(3);

    // Solve for X
    X_temp = A1.inverse() * (B1 - C1 * Y_temp - C2 * Z_temp);
    X_final = X_temp;

    // Solve for Y
    Eigen::Matrix4d A1p = A1 - C1 * A2.inverse() * C1.transpose() * opt * nstd1;
    Eigen::Matrix4d B1p = B1 - C1 * A2.inverse() * B2 * opt * nstd1;
    Y_temp = A1p.inverse() * B1p;
    Y_final = Y_temp;

    // Solve for Z
    Eigen::Matrix4d A2p = A2 - C2 * A1.inverse() * C2.transpose() * opt * nstd2;
    Eigen::Matrix4d B2p = B2 - C2 * A1.inverse() * B1 * opt * nstd2;
    Z_temp = A2p.inverse() * B2p;
    Z_final = Z_temp;
}

