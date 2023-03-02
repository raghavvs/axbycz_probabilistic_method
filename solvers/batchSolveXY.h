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

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

Matrix4d vec2mat(VectorXd vec) {
    Matrix4d mat;
    mat << vec(0), vec(1), vec(2), vec(3),
           vec(4), vec(5), vec(6), vec(7),
           vec(8), vec(9), vec(10), vec(11),
           vec(12), vec(13), vec(14), vec(15);
    return mat;
}

VectorXd mat2vec(Matrix4d mat) {
    VectorXd vec(16);
    vec << mat(0,0), mat(0,1), mat(0,2), mat(0,3),
           mat(1,0), mat(1,1), mat(1,2), mat(1,3),
           mat(2,0), mat(2,1), mat(2,2), mat(2,3),
           mat(3,0), mat(3,1), mat(3,2), mat(3,3);
    return vec;
}

void distibutionProps(const MatrixXd& data, VectorXd& mean, MatrixXd& cov) {
    int num_samples = data.cols();
    mean = data.rowwise().mean();
    MatrixXd centered = data.colwise() - mean;
    cov = centered * centered.transpose() / (num_samples - 1);
}

void batchSolveXY(Matrix4d A[8], Matrix4d B[8], bool opt, double nstd_A, double nstd_B, Matrix4d X_candidate[8], Matrix4d Y_candidate[8]) {
    MatrixXd A_mex(4, 8);
    MatrixXd B_mex(4, 8);

    for (int i = 0; i < 8; i++) {
        A_mex.col(i) = mat2vec(A[i]);
        B_mex.col(i) = mat2vec(B[i]);
    }

    VectorXd MeanA, MeanB;
    MatrixXd SigA, SigB;

    distibutionProps(A_mex, MeanA, SigA);
    distibutionProps(B_mex, MeanB, SigB);

    if (opt) {
        SigA -= nstd_A * MatrixXd::Identity(6, 6);
        SigB -= nstd_B * MatrixXd::Identity(6, 6);
    }

    SelfAdjointEigenSolver<MatrixXd> eigA(SigA.block<3,3>(0,0));
    SelfAdjointEigenSolver<MatrixXd> eigB(SigB.block<3,3>(0,0));
    MatrixXd VA = eigA.eigenvectors();
    MatrixXd VB = eigB.eigenvectors();

    MatrixXd Q1 = MatrixXd::Identity(3, 3);
    MatrixXd Q2(3, 3), Q3(3, 3), Q4(3, 3);
    Q2 << -1, 0, 0,
           0, -1, 0,
           0, 0, 1;
    Q3 << -1, 0, 0,
           0, 1, 0,
           0, 0, -1;
    Q4 << 1, 0, 0,
           0, -1, 0,
           0, 0, -1;

    MatrixXd Rx_solved[8];
    for (int i = 0; i < 8; i++) {
        MatrixXd tx_temp = ((Rx_solved[i].transpose()*SigA.block<3,3>(0,0)*Rx_solved[i]).inverse() *
                            (SigB.block<3,3>(0,3)-Rx_solved[i].transpose()*SigA.block<3,3>(0,3)*Rx_solved[i])).transpose();
        MatrixXd tx = -Rx_solved[i] * tx_temp;
        X_candidate[i] << Rx_solved[i], tx, 0, 0, 0, 1;
        Y_candidate[i] = MeanA * X_candidate[i] / (MeanB);
    }
    X = X_candidate;
    Y = Y_candidate;
}