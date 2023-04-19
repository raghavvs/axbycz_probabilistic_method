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
#include <Eigen/Eigenvalues>
#include <numeric>
#include "meanCov.h"
#include "so3Vec.h"

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

    std::vector<Eigen::Matrix4d> X_candidate(8, Eigen::Matrix4d::Zero());
    std::vector<Eigen::Matrix4d> Y_candidate(8, Eigen::Matrix4d::Zero());

    // Calculate mean and covariance for A and B
    meanCov(A, MeanA, SigA);
    meanCov(B, MeanB, SigB);

    // update SigA and SigB if nstd_A and nstd_B are known
    if (opt) {
        SigA -= nstd_A * Eigen::MatrixXd::Identity(6, 6);
        SigB -= nstd_B * Eigen::MatrixXd::Identity(6, 6);
    }

    /*// Calculate eigenvectors of top left 3x3 sub-matrices of SigA and SigB
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_solver_A(SigA.topLeftCorner<3, 3>());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_solver_B(SigB.topLeftCorner<3, 3>());*/

    /*Eigen::MatrixXd A_temp = SigA.block<3, 3>(0, 0);
    Eigen::MatrixXd B_temp = SigB.block<3, 3>(0, 0);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_A(A_temp, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_B(B_temp, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::MatrixXd VA = svd_A.matrixU() * svd_A.matrixV().transpose();
    Eigen::MatrixXd VB = svd_B.matrixU() * svd_B.matrixV().transpose();*/

    /*Eigen::MatrixXd A_temp = SigA.block<3, 3>(0, 0);
    Eigen::MatrixXd B_temp = SigB.block<3, 3>(0, 0);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver_A(A_temp);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver_B(B_temp);

    Eigen::MatrixXd VA = eig_solver_A.eigenvectors();
    Eigen::MatrixXd VB = eig_solver_B.eigenvectors();*/

    Eigen::EigenSolver<Eigen::MatrixXd> esA(SigA.block<3, 3>(0, 0));
    Eigen::MatrixXd VA = esA.eigenvectors().real();

    Eigen::EigenSolver<Eigen::MatrixXd> esB(SigB.block<3, 3>(0, 0));
    Eigen::MatrixXd VB = esB.eigenvectors().real();

    std::cout << "Eigenvectors A (VA): " << std::endl;
    std::cout << VA << std::endl;
    std::cout << "Eigenvectors B (VB): " << std::endl;
    std::cout << VB << std::endl;

    // Define Q matrices
    Eigen::Matrix3d Q1, Q2, Q3, Q4;
    Q1 = Eigen::Matrix3d::Identity();
    Q2 = (Eigen::Matrix3d() << -1, 0, 0, 0, -1, 0, 0, 0, 1).finished();
    Q3 = (Eigen::Matrix3d() << -1, 0, 0, 0, 1, 0, 0, 0, -1).finished();
    Q4 = (Eigen::Matrix3d() << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished();

    std::vector<Eigen::Matrix3d> Rx_solved(8, Eigen::Matrix3d::Zero());

    // There are eight possibilities for Rx
    Rx_solved[0] = VA * Q1 * VB.transpose();
    Rx_solved[1] = VA * Q2 * VB.transpose();
    Rx_solved[2] = VA * Q3 * VB.transpose();
    Rx_solved[3] = VA * Q4 * VB.transpose();
    Rx_solved[4] = VA * (-Q1) * VB.transpose();
    Rx_solved[5] = VA * (-Q2) * VB.transpose();
    Rx_solved[6] = VA * (-Q3) * VB.transpose();
    Rx_solved[7] = VA * (-Q4) * VB.transpose();

    std::cout << "Rx_solved[0]: " << std::endl;
    std::cout << Rx_solved[0] << std::endl;
    std::cout << "Rx_solved[4]: " << std::endl;
    std::cout << Rx_solved[4] << std::endl;

    X.resize(8);
    Y.resize(8);

    for (int i = 0; i < 8; ++i) {
        Eigen::Matrix3d temp = (Rx_solved[i].transpose() * SigA.block<3, 3>(0, 0) * Rx_solved[i]).inverse() *
                               (SigB.block<3, 3>(0, 3) - Rx_solved[i].transpose() * SigA.block<3, 3>(0, 3) * Rx_solved[i]);

        Eigen::Vector3d tx_temp = so3Vec(temp.transpose());
        Eigen::Vector3d tx = -Rx_solved[i] * tx_temp;

        X_candidate[i] << Rx_solved[i], tx, Eigen::Vector3d::Zero().transpose(), 1;
        Y_candidate[i] = MeanA * X_candidate[i] * MeanB.inverse();

        // Set the output X and Y
        X[i] = X_candidate[i];
        Y[i] = Y_candidate[i];
    }
}

int main() {
    // Create deterministic input matrices A and B
    std::vector<Eigen::Matrix4d> A(2), B(2);

    A[0] << 0.9363, -0.2751, 0.2183, 1.2020,
            0.2896, 0.9566, -0.0392, -0.1022,
            -0.1985, 0.0978, 0.9750, 0.3426,
            0.0, 0.0, 0.0, 1.0;

    A[1] << 0.9938, -0.0975, 0.0599, -0.2246,
            0.0975, 0.9951, -0.0273, 0.1088,
            -0.0603, 0.0250, 0.9981, 0.4839,
            0.0, 0.0, 0.0, 1.0;

    B[0] << 0.8660, -0.2896, 0.4082, 0.9501,
            0.5000, 0.8660, -0.0000, -0.5507,
            -0.0000, 0.0000, 1.0000, 0.5000,
            0.0, 0.0, 0.0, 1.0;

    B[1] << 0.9603, -0.1944, 0.2014, 0.6231,
            0.2791, 0.6829, -0.6752, -0.4567,
            -0.0000, 0.7071, 0.7071, 0.7071,
            0.0, 0.0, 0.0, 1.0;

    bool opt = false;
    double nstd_A = 0.0, nstd_B = 0.0;
    std::vector<Eigen::Matrix4d> X, Y;
    Eigen::Matrix4d MeanA, MeanB;
    Eigen::Matrix<double, 6, 6> SigA, SigB;

    batchSolveXY(A, B, opt, nstd_A, nstd_B, X, Y, MeanA, MeanB, SigA, SigB);

    // Display results
    std::cout << "X:\n";
    for (const auto &x: X) {
        std::cout << x << "\n\n";
    }

    std::cout << "Y:\n";
    for (const auto &y: Y) {
        std::cout << y << "\n\n";
    }

    std::cout << "MeanA:\n" << MeanA << "\n\n";
    std::cout << "MeanB:\n" << MeanB << "\n\n";
    std::cout << "SigA:\n" << SigA << "\n\n";
    std::cout << "SigB:\n" << SigB << "\n\n";

    return 0;
}
