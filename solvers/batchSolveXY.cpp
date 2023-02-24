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
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <so3Vec.h>

using namespace std;
using namespace Eigen;

Matrix<double, 4, 4, RowMajor> distibutionPropsMex(const MatrixXd &M)
{
    int r = M.rows();
    int c = M.cols();

    MatrixXd M_centered = M.colwise() - M.rowwise().mean();
    MatrixXd S = (M_centered * M_centered.transpose()) / double(c - 1);

    Matrix<double, 4, 4, RowMajor> M_props;
    M_props.topLeftCorner<3, 3>() = S.block<3, 3>(0, 0);
    M_props.topRightCorner<3, 1>() = S.block<0, 3>(3, 0);
    M_props.bottomLeftCorner<1, 3>().setZero();
    M_props.bottomRightCorner<1, 1>().setZero();
    M_props(3, 3) = 1.0;

    return M_props;
}

void batchSolveXY(MatrixXd &A, MatrixXd &B, bool opt, double nstd_A, double nstd_B,
    Matrix<double, 4, 4, RowMajor> &MeanA, Matrix<double, 4, 4, RowMajor> &MeanB,
    Matrix<double, 6, 6, RowMajor> &SigA, Matrix<double, 6, 6, RowMajor> &SigB,
    Matrix<double, 4, 4, RowMajor> *X, Matrix<double, 4, 4, RowMajor> *Y)
{
    Matrix<double, 4, 4, RowMajor> X_candidate[8];
    Matrix<double, 4, 4, RowMajor> Y_candidate[8];

    // Reshape A and B for matching the input sizes of mex functions
    int a1 = A.rows();
    int a2 = A.cols();
    int a3 = A.size() / (a1 * a2);

    MatrixXd SigA_13 = SigA.block(0, 0, 3, 3);
    SigA_13.resize(1, 9);
    MatrixXd SigB_13 = SigA.block(0, 0, 3, 3);
    SigB_13.resize(1, 9);

    // [MeanA, SigA] = distibutionPropsMex(A_mex);
    // [MeanB, SigB] = distibutionPropsMex(B_mex);
    MeanA = distibutionPropsMex(A_mex);
    MeanB = distibutionPropsMex(B_mex);
    SigA.setZero();
    SigB.setZero();
    for (int i = 0; i < a3; i++)
    {
        MatrixXd A_i = A.block<4, 4>(0, 0, a1, a2).transpose();
        MatrixXd B_i = B.block<4, 4>(0, 0, a1, a2).transpose();
        MatrixXd A_i_centered = A_i.colwise() - MeanA.transpose();
        MatrixXd B_i_centered = B_i.colwise() - MeanB.transpose();
        SigA += (A_i_centered * A_i_centered.transpose()) / double(a3 - 1);
        SigB += (B_i_centered * B_i_centered.transpose()) / double(a3 - 1);
    }
    SigA /= double(a3);
    SigB /= double(a3);

    // update SigA and SigB if nstd_A and nstd_B are known
    if (opt) {
        SigA -= nstd_A * MatrixXd::Identity(6, 6);
        SigB -= nstd_B * MatrixXd::Identity(6, 6);
        }

    MatrixXd VA(3, 3), VB(3, 3);
    MatrixXd SigA_13 = SigA.block<3, 3>(0, 0);
    MatrixXd SigB_13 = SigB.block<3, 3>(0, 0);

    SelfAdjointEigenSolver<MatrixXd> eigA(SigA_13), eigB(SigB_13);
    MatrixXd VA = eigA.eigenvectors();
    MatrixXd VB = eigB.eigenvectors();

    Matrix3d Q1, Q2, Q3, Q4;
    Q1 = Matrix3d::Identity();
    Q2 << -1, 0, 0,
    0,-1, 0,
    0, 0, 1;
    Q3 << -1, 0, 0,
    0, 1, 0,
    0, 0,-1;
    Q4 << 1, 0, 0,
    0,-1, 0,
    0, 0,-1;

    Matrix3d Rx_solved[8];
    Rx_solved[0] = VA * Q1 * VB.transpose();
    Rx_solved[1] = VA * Q2 * VB.transpose();
    Rx_solved[2] = VA * Q3 * VB.transpose();
    Rx_solved[3] = VA * Q4 * VB.transpose();
    Rx_solved[4] = VA * -Q1 * VB.transpose();
    Rx_solved[5] = VA * -Q2 * VB.transpose();
    Rx_solved[6] = VA * -Q3 * VB.transpose();
    Rx_solved[7] = VA * -Q4 * VB.transpose();

    MatrixXd SigA_13_inv = (Rx_solved[0].transpose() * SigA_13 * Rx_solved[0]).inverse();
    MatrixXd SigB_14_16 = SigB.block(0, 3, 3, 3) - Rx_solved[0].transpose() * SigA.block(0, 3, 3, 3) * Rx_solved[0];
    MatrixXd tx_temp = -SigA_13_inv * SigB_14_16;
    Vector3d tx = -Rx_solved[0] * so3Vec(tx_temp).transpose();

    Matrix4d X_candidate[8], Y_candidate[8];
    for (int i = 0; i < 8; i++) {
        X_candidate[i].block(0, 0, 3, 3) = Rx_solved[i];
        X_candidate[i].block(0, 3, 3, 1) = tx;
        X_candidate[i].row(3) << 0, 0, 0, 1;
        Y_candidate[i] = MeanA * X_candidate[i] * MeanB.inverse();
    }

    for (int i = 0; i < 8; i++) {
        X[i] = X_candidate[i];
        Y[i] = Y_candidate[i];
    }
}