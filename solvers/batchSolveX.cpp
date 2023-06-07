/*
Refer Quanli Ma's GitHub repo for original MATLAB code
 AX = XB - Solver
*/

#include <iostream>
#include <algorithm>
#include <vector>
#include <Eigen/Eigenvalues>
#include "meanCov.h"
#include "so3Vec.h"
#include "paramExtract.h"

// Sorting function
void sortEigenVectors(const Eigen::VectorXd& eigenvalues,
                      const Eigen::MatrixXd& eigenvectors,
                      Eigen::VectorXd& sorted_eigenvalues,
                      Eigen::MatrixXd& sorted_eigenvectors) {
    // Create a vector of pairs (eigenvalue, eigenvector)
    std::vector<std::pair<double, Eigen::VectorXd>> eigen_pairs(eigenvalues.size());

    for (size_t i = 0; i < eigenvalues.size(); ++i) {
        eigen_pairs[i] = std::make_pair(eigenvalues[i], eigenvectors.col(i));
    }

    // Sort the eigen_pairs based on eigenvalues
    std::sort(eigen_pairs.begin(), eigen_pairs.end(), [](const std::pair<double, Eigen::VectorXd>& a, const std::pair<double, Eigen::VectorXd>& b) {
        return a.first < b.first;
    });

    // Fill the sorted_eigenvalues and sorted_eigenvectors with the sorted data
    for (size_t i = 0; i < eigen_pairs.size(); ++i) {
        sorted_eigenvalues[i] = eigen_pairs[i].first;
        sorted_eigenvectors.col(i) = eigen_pairs[i].second;
    }
}

void batchSolveXY(const std::vector<Eigen::Matrix4d> &A,
                  const std::vector<Eigen::Matrix4d> &B,
                  bool opt,
                  std::vector<Eigen::Matrix4d> &X,
                  Eigen::Matrix4d &MeanA,
                  Eigen::Matrix4d &MeanB,
                  Eigen::Matrix<double, 6, 6> &SigA,
                  Eigen::Matrix<double, 6, 6> &SigB,
                  double t_error) {

    std::vector<Eigen::Matrix4d> X_candidate(8, Eigen::Matrix4d::Zero());

    // Calculate mean and covariance for A and B
    meanCov(A, MeanA, SigA);
    meanCov(B, MeanB, SigB);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esA(SigA.block<3, 3>(0, 0));
    Eigen::MatrixXd VA = esA.eigenvectors().real();
    Eigen::VectorXcd eigenvalues_A = esA.eigenvalues(); // Eigenvalues for SigA block

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esB(SigB.block<3, 3>(0, 0));
    Eigen::MatrixXd VB = esB.eigenvectors().real();
    Eigen::VectorXcd eigenvalues_B = esB.eigenvalues(); // Eigenvalues for SigB block

    Eigen::VectorXd real_eigenvalues_A = eigenvalues_A.real();
    Eigen::VectorXd real_eigenvalues_B = eigenvalues_B.real();

    // Sort eigenvalues and eigenvectors for SigA block
    Eigen::MatrixXd sorted_VA(3, 3);
    Eigen::VectorXd sorted_eigenvalues_A(3);
    sortEigenVectors(real_eigenvalues_A, VA, sorted_eigenvalues_A, sorted_VA);

    // Sort eigenvalues and eigenvectors for SigB block
    Eigen::MatrixXd sorted_VB(3, 3);
    Eigen::VectorXd sorted_eigenvalues_B(3);
    sortEigenVectors(real_eigenvalues_B, VB, sorted_eigenvalues_B, sorted_VB);

    // Define Q matrices
    Eigen::Matrix3d Q1, Q2, Q3, Q4;
    Q1 = Eigen::Matrix3d::Identity();
    Q2 = (Eigen::Matrix3d() << -1, 0, 0, 0, -1, 0, 0, 0, 1).finished();
    Q3 = (Eigen::Matrix3d() << -1, 0, 0, 0, 1, 0, 0, 0, -1).finished();
    Q4 = (Eigen::Matrix3d() << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished();

    std::vector<Eigen::Matrix3d> Rx_solved(8, Eigen::Matrix3d::Zero());

    // There are eight possibilities for Rx
    Rx_solved[0] = sorted_VA * Q1 * sorted_VB.transpose();
    Rx_solved[1] = sorted_VA * Q2 * sorted_VB.transpose();
    Rx_solved[2] = sorted_VA * Q3 * sorted_VB.transpose();
    Rx_solved[3] = sorted_VA * Q4 * sorted_VB.transpose();
    Rx_solved[4] = sorted_VA * (-Q1) * sorted_VB.transpose();
    Rx_solved[5] = sorted_VA * (-Q2) * sorted_VB.transpose();
    Rx_solved[6] = sorted_VA * (-Q3) * sorted_VB.transpose();
    Rx_solved[7] = sorted_VA * (-Q4) * sorted_VB.transpose();

    double ta, tb, da, db;
    Eigen::MatrixXd Na, Nb,Xa, Xb;
    paramExtract(ta, Na, da, pa, Xa);
    paramExtract(tb, Nb, db, pb, Xb);

    Eigen::VectorXd na = so3Vec(Na);
    Eigen::VectorXd nb = so3Vec(Nb);

    X.resize(8);

    double min = std::numeric_limits<double>::infinity();
    Eigen::MatrixXd Rx;

    // Iterate over Rx_solved
    for (int i = 0; i < 8; ++i) {
        Eigen::MatrixXd curr_Rx = Rx_solved.block<3,3>(0,0,i);  // Assuming Rx_solved is a 3D matrix
        if (std::abs(curr_Rx.determinant() - 1) < 0.001 && (na - curr_Rx * nb).norm() < min) {
            min = (na - curr_Rx * nb).norm();
            Rx = curr_Rx;
        }
    }

    // Compute tx
    Eigen::VectorXd tx_temp = so3_vec(((Rx.transpose() * SigA.block<3,3>(0,0) * Rx).inverse() * (SigB.block<3,3>(0,0) - Rx.transpose() * SigA.block<3,3>(0,0) * Rx)).transpose());
    Eigen::VectorXd tx = -Rx * tx_temp;

    // Compute X
    Eigen::MatrixXd X(4, 4);
    X << Rx, tx,
            0, 0, 0, 1;

    // Compute t_error
    Eigen::VectorXd t_error = (MeanA.block<3,3>(0,0) - Eigen::MatrixXd::Identity(3,3)) * tx - Rx * MeanB.block<3,1>(0,4) + MeanA.block<3,1>(0,4);
    double error_norm = t_error.norm();
}