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

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "meanCov.h"
#include "so3Vec.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>

void batchSolveXY(const std::vector<Eigen::Matrix4d> &A,
                  const std::vector<Eigen::Matrix4d> &B,
                  bool opt,
                  double nstd_A,
                  double nstd_B,
                  std::vector<Eigen::Matrix4d> &X,
                  std::vector<Eigen::Matrix4d> &Y,
                  Eigen::Matrix4d &MeanA,
                  Eigen::Matrix4d &MeanB,
                  Eigen::Matrix<double, 6, 6> &SigA,
                  Eigen::Matrix<double, 6, 6> &SigB) {

    X.resize(8);
    Y.resize(8);

    // Compute Mean and Covariance for A and B
    meanCov(A, MeanA, SigA);
    meanCov(B, MeanB, SigB);

    // Update SigA and SigB if nstd_A and nstd_B are known
    if (opt) {
        SigA -= nstd_A * Eigen::MatrixXd::Identity(6, 6);
        SigB -= nstd_B * Eigen::MatrixXd::Identity(6, 6);
    }

    // Compute eigenvectors of SigA and SigB
    Eigen::Matrix3d VA = SigA.block<3, 3>(0, 0).eigenvalues().real().asDiagonal();
    Eigen::Matrix3d VB = SigB.block<3, 3>(0, 0).eigenvalues().real().asDiagonal();

    // Define the Q matrices
    Eigen::Matrix3d Q1 = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d Q2; Q2 << -1, 0, 0, 0, -1, 0, 0, 0, 1;
    Eigen::Matrix3d Q3; Q3 << -1, 0, 0, 0, 1, 0, 0, 0, -1;
    Eigen::Matrix3d Q4; Q4 << 1, 0, 0, 0, -1, 0, 0, 0, -1;

    std::vector<Eigen::Matrix3d> Rx_solved(8);

    // Compute 8 possible Rx matrices
    for (int i = 0; i < 8; i++) {
        Rx_solved[i] = VA * (i < 4 ? Q1 : -Q1) * VB.transpose();
        Q1 = (i % 4 == 0) ? Q2 : ((i % 4 == 1) ? Q3 : ((i % 4 == 2) ? Q4 : Q1));
    }

    // Compute X and Y candidates
    for (int i = 0; i < 8; i++) {
        Eigen::Matrix3d temp = ((Rx_solved[i].transpose() * SigA.block<3, 3>(0, 0) *
                                Rx_solved[i]).inverse() * (SigB.block<3, 3>(0, 3) -
                                Rx_solved[i].transpose() * SigA.block<3, 3>(0, 3) *
                                Rx_solved[i])).transpose();
        Eigen::Vector3d tx_temp = so3Vec(temp.transpose());
        Eigen::Vector3d tx = -Rx_solved[i] * tx_temp;
        X[i] << Rx_solved[i], tx, 0, 0, 0, 1;
        Y[i] = MeanA * X[i] * MeanB.inverse();
    }
}


int main() {
    int len = 10;
    bool opt = true;
    double nstd_A = 0.1;
    double nstd_B = 0.2;
    srand(12345);

    std::vector<Eigen::Matrix4d> A(len), B(len);
    for(int i = 0; i < len; i++){
        A[i] = Eigen::Matrix4d::Random();
        B[i] = Eigen::Matrix4d::Random();
    }

    std::vector<Eigen::Matrix4d> X(len), Y(len);
    Eigen::Matrix4d MeanA, MeanB;
    Eigen::Matrix<double, 6, 6> SigA, SigB;

    batchSolveXY(A, B, opt, nstd_A, nstd_B,
                 X,Y,
                 MeanA,
                 MeanB,
                 SigA,
                 SigB);

    std::cout << "X: \n" << X[0] << std::endl;
    std::cout << "Y: \n" << Y[0] << std::endl;

    return 0;
}