/*
The code defines two functions: distributionPropsMex and batchSolveXY.
The former function takes a matrix as an input and calculates the
covariance matrix of the matrix along with other properties, which
it stores in a matrix that it returns. The batchSolveXY function
takes several inputs and solves for the rotation matrix X and Y.
It does this by first calling the distributionPropsMex function for
two matrices A and B, to obtain their mean and covariance matrix.
It then calculates the eigenvectors of the covariance matrices,
and uses them to calculate eight possible solutions for the
rotation matrix. It stores these solutions in two arrays, one
for X and one for Y.

The batchSolveXY function can also adjust the covariance matrices
of A and B based on nstd_A and nstd_B. If the boolean input "opt"
is true, then the covariance matrices are adjusted by subtracting
the identity matrix multiplied by nstd_A and nstd_B. The code does
not explain what these values are or what they represent, so it
is unclear what effect this has on the calculation. Additionally,
the code contains some errors, such as redefining SigA_13 and
SigB_13, which causes a compiler error, and using an ellipsis
(...) instead of an integer to index into the Rx_solved array.

Input:
    A, B: Matrices - dim - 4x4 - pass by reference
    len: number of data pairs / A, B matrices
    opt: update SigA and SigB if nstd_A and nstd_B are known
    nstd_A, nstd_B: standard deviation of A, B matrices
Output:
    X, Y: Matrices - dim - 4x4
    MeanA, MeanB: Matrices - dim - 4x4 - Mean of A, B
    SigA, SigB: Matrices - dim - 6x6 - Covariance of A, B
*/

#ifndef BATCHSOLVEXY_H
#define BATCHSOLVEXY_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "meanCov.h"
#include "so3Vec.h"

void batchSolveXY(const Eigen::Matrix4d& A,
                  const Eigen::Matrix4d& B,
                  int len,
                  bool opt,
                  double nstd_A,
                  double nstd_B,
                  std::vector<Eigen::MatrixXd> &X,
                  std::vector<Eigen::MatrixXd> &Y,
                  Eigen::MatrixXd& MeanA,
                  Eigen::MatrixXd& MeanB,
                  Eigen::MatrixXd& SigA,
                  Eigen::MatrixXd& SigB) {

    std::vector<Eigen::Matrix4d> X_candidate(8), Y_candidate(8);

    // Calculate mean and covariance for A and B
    meanCov(A, len, MeanA, SigA);
    meanCov(B, len, MeanB, SigB);

    // update SigA and SigB if nstd_A and nstd_B are known
    if (opt) {
        SigA -= nstd_A * Eigen::MatrixXd::Identity(6, 6);
        SigB -= nstd_B * Eigen::MatrixXd::Identity(6, 6);
    }

    // Calculate eigenvectors of top left 3x3 sub-matrices of SigA and SigB
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_solver_A(SigA.topLeftCorner<3, 3>());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_solver_B(SigB.topLeftCorner<3, 3>());

    auto const& VA = eig_solver_A.eigenvectors();
    auto const& VB = eig_solver_B.eigenvectors();

    // Define Q matrices
    Eigen::MatrixXd Q1, Q2, Q3, Q4;
    Q1 = Eigen::MatrixXd::Identity(3, 3);
    Q2 = (Eigen::MatrixXd(3, 3) << -1, 0, 0, 0, -1, 0, 0, 0, 1).finished();
    Q3 = (Eigen::MatrixXd(3, 3) << -1, 0, 0, 0, 1, 0, 0, 0, -1).finished();
    Q4 = (Eigen::MatrixXd(3, 3) << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished();

    Eigen::Matrix3d Rx_solved[8];

    // There are eight possibilities for Rx
    Rx_solved[0] = VA * Q1 * VB.transpose();
    Rx_solved[1] = VA * Q2 * VB.transpose();
    Rx_solved[2] = VA * Q3 * VB.transpose();
    Rx_solved[3] = VA * Q4 * VB.transpose();
    Rx_solved[4] = VA * (-Q1) * VB.transpose();
    Rx_solved[5] = VA * (-Q2) * VB.transpose();
    Rx_solved[6] = VA * (-Q3) * VB.transpose();
    Rx_solved[7] = VA * (-Q4) * VB.transpose();

    X.resize(8);
    Y.resize(8);

    for (int i = 0; i < 8; i++) {
        // block SigA and SigB to 3x3 sub-matrices
        Eigen::Matrix3d sigA_33 = SigA.block<3, 3>(0, 0);
        Eigen::Matrix3d sigB_33 = SigB.block<3, 3>(0, 3);

        Eigen::Matrix3d temp = (Rx_solved[i].transpose() * sigA_33 * Rx_solved[i]).inverse() *
                               (sigB_33 - Rx_solved[i].transpose() * SigA.block<3, 3>(0, 3) * Rx_solved[i]);

        Eigen::Vector3d tx_temp = so3Vec(temp.transpose());

        // Construct X and Y candidates
        X_candidate[i] << Rx_solved[i], -Rx_solved[i] * tx_temp, Eigen::Vector4d::Zero().transpose();
        Y_candidate[i] = MeanA * X_candidate[i] * MeanB.inverse();

        // Set the output X and Y
        X[i] = X_candidate[i];
        Y[i] = Y_candidate[i];
    }

    // Set the output MeanA, MeanB, SigA, and SigB
    for (int i = 0; i < 8; i++) {
        MeanA = MeanA * X[i] * MeanB.inverse();
    }
    MeanB = Eigen::Matrix4d::Identity();
    SigA = SigA.block<3, 3>(0, 0);
    SigB = SigB.block<3, 3>(0, 3);
}

#endif