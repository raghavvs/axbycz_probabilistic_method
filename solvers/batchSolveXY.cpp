/*
The code defines two functions: distibutionPropsMex and batchSolveXY. 
The former function takes a matrix as an input and calculates the 
covariance matrix of the matrix along with other properties, which 
it stores in a matrix that it returns. The batchSolveXY function 
takes several inputs and solves for the rotation matrix X and Y. 
It does this by first calling the distibutionPropsMex function for 
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
#include <Eigen/Dense>
#include <vector>
#include <meanCov.h>

using namespace Eigen;

std::vector<MatrixXd> eig(MatrixXd A) {
    EigenSolver<MatrixXd> es(A);
    std::vector<MatrixXd> result(2);
    result[0] = es.eigenvalues().real();
    result[1] = es.eigenvectors().real();
    return result;
}

void batchSolveXY(MatrixXd A, MatrixXd B, bool opt, double nstd_A, double nstd_B, MatrixXd& X, MatrixXd& Y) {
    X = MatrixXd::Zero(4, 4, 8);
    Y = MatrixXd::Zero(4, 4, 8);

    int a1 = A.rows(), a2 = A.cols(), a3 = A.size() / (a1 * a2);
    MatrixXd A_mex = A.reshape(a1, a2 * a3);
    MatrixXd B_mex = B.reshape(a1, a2 * a3);

    MatrixXd MeanA, SigA, MeanB, SigB;
    MeanA.setZero(4, 1);
    SigA.setZero(6, 6);
    MeanB.setZero(4, 1);
    SigB.setZero(6, 6);

    MeanA = meanCov(A_mex).first;
    SigA = meanCov(A_mex).second;
    MeanB = meanCov(B_mex).first;
    SigB = meanCov(B_mex).second;

    if (opt) {
        SigA -= nstd_A * MatrixXd::Identity(6, 6);
        SigB -= nstd_B * MatrixXd::Identity(6, 6);
    }

    MatrixXd VA, VB;
    VA.setZero(3, 3);
    VB.setZero(3, 3);

    std::vector<MatrixXd> VA_eig, VB_eig;
    VA_eig = eig(SigA.block<3, 3>(0, 0));
    VB_eig = eig(SigB.block<3, 3>(0, 0));
    VA = VA_eig[1];
    VB = VB_eig[1];

    MatrixXd Q1, Q2, Q3, Q4;
    Q1 = MatrixXd::Identity(3, 3);
    Q2 = (MatrixXd(3, 3) << -1, 0, 0, 0, -1, 0, 0, 0, 1).finished();
    Q3 = (MatrixXd(3, 3) << -1, 0, 0, 0, 1, 0, 0, 0, -1).finished();
    Q4 = (MatrixXd(3, 3) << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished();

    MatrixXd Rx_solved(3, 3, 8);

    // There are 8 possibilities of Rx
    Rx_solved.slice(0) = VA * Q1 * VB.transpose();
    Rx_solved.slice(1) = VA * Q2 * VB.transpose();
    Rx_solved.slice(2) = VA * Q3 * VB.transpose();
    Rx_solved.slice(3) = VA * Q4 * VB.transpose();
    Rx_solved.slice(4) = VA * -Q1 * VB.transpose();
    Rx_solved.slice(5) = VA * -Q2 * VB.transpose();
    Rx_solved.slice(6) = VA * -Q3 * VB.transpose();
    Rx_solved.slice(7) = VA * -Q4 * VB.transpose();

    for (int i = 0; i < 8; i++) {
        MatrixXd SigA_sub = SigA.block(0, 0, 3, 3);
        MatrixXd SigB_sub = SigB.block(0, 0, 3, 3);

        VectorXd tx_temp_vec = (Rx_solved.slice(i).transpose() * SigA_sub * Rx_solved.slice(i)).inverse() * 
                            (SigB_sub - Rx_solved.slice(i).transpose() * SigA.block(0, 3, 3, 3));
        Vector3d tx_temp(tx_temp_vec(0), tx_temp_vec(1), tx_temp_vec(2));
        Vector3d tx = -Rx_solved.slice(i) * tx_temp;

        Matrix4d X_candidate;
        X_candidate.block(0, 0, 3, 3) = Rx_solved.slice(i);
        X_candidate.block(0, 3, 3, 1) = tx;
        X_candidate.block(3, 0, 1, 3) = MatrixXd::Zero(1, 3);
        X_candidate(3, 3) = 1;

        Matrix4d Y_candidate = MeanA * X_candidate / MeanB;

        X.slice(i) = X_candidate;
        Y.slice(i) = Y_candidate;
    }
}