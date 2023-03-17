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
*/

#include <iostream>
#include <array>
#include <Eigen/Dense>
#include "meanCov.h"
#include "so3Vec.h"
#include "se3Vec.h"

void batchSolveXY(const Eigen::Matrix4d* A,
                  const Eigen::Matrix4d* B,
                  int len,
                  bool opt,
                  double nstd_A,
                  double nstd_B,
                  Eigen::Matrix4d& X,
                  Eigen::Matrix4d& Y,
                  Eigen::Matrix4d& MeanA,
                  Eigen::Matrix4d& MeanB,
                  Eigen::MatrixXd& SigA,
                  Eigen::MatrixXd& SigB) {

    Eigen::Matrix4d X_candidate, Y_candidate;

    // Reshape A and B for matching the input sizes of mex functions
    int a1 = A[0].rows();
    int a2 = A[0].cols();
    Eigen::Matrix3d A_mex = Eigen::Map<Eigen::Matrix3d>(A->data(), a1, a2 * len);
    Eigen::Matrix3d B_mex = Eigen::Map<Eigen::Matrix3d>(B->data(), a1, a2 * len);

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

    auto VA = eig_solver_A.eigenvectors();
    auto VB = eig_solver_B.eigenvectors();

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

    for (int i = 0; i < 8; i++) {
        Eigen::MatrixXd SigA_sub = SigA.block(0, 0, 3, 3);
        Eigen::MatrixXd SigB_sub = SigB.block(0, 0, 3, 3);

        Eigen::VectorXd tx_temp_vec = (Rx_solved.block<3, 3>(0, i).transpose() * SigA_sub * Rx_solved.block<3, 3>(0, i)).inverse() *
                               (SigB_sub - Rx_solved.block<3, 3>(0, i).transpose() * SigA.block(0, 3, 3, 3));
        Eigen::Vector3d tx_temp(tx_temp_vec(0), tx_temp_vec(1), tx_temp_vec(2));
        Eigen::Vector3d tx = -Rx_solved.block<3, 3>(0, i) * tx_temp;

        X_candidate.block(0, 0, 3, 3) = Rx_solved.block<3, 3>(0, i);
        X_candidate.block(0, 3, 3, 1) = tx;
        X_candidate.block(3, 0, 1, 3) = Eigen::MatrixXd::Zero(1, 3);
        X_candidate(3, 3) = 1;

        Eigen::Matrix4d Y_candidate = MeanA * X_candidate * MeanB.inverse();

        X.block(0, 0, 4, 4) = X_candidate;
        Y.block(0, 0, 4, 4) = Y_candidate;
    }
    X = X_candidate;
    Y = Y_candidate;
}


/*
for (int i = 0; i < 8; i++) {
Eigen::Matrix3d temp = (Rx_solved[i].transpose() * SigA.topLeftCorner<3, 3>() * Rx_solved[i]).inverse() *
                       (SigB.topRightCorner<3, 3>() - Rx_solved[i].transpose() * SigA.topRightCorner<3, 3>() * Rx_solved[i]);
Eigen::Vector3d tx_temp = se3Vec(temp);
Eigen::Vector3d tx = -Rx_solved[i] * tx_temp;
X[i].block<3,3>(0,0) = Rx_solved[i];
X[i].block<3,1>(0,3) = tx;
X[i].block<1,4>(3,0) << 0, 0 , 0 ,1;
Y[i] = MeanA * X[i] / MeanB;
}*/
