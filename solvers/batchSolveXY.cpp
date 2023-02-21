#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

void batchSolveXY(Matrix<float, 4, 4, Eigen::RowMajor> *A, Matrix<float, 4, 4, Eigen::RowMajor> *B, bool opt, float nstd_A, float nstd_B,
                  Matrix<float, 4, 4, Eigen::RowMajor> *X, Matrix<float, 4, 4, Eigen::RowMajor> *Y,
                  Matrix<float, 3, 1> *MeanA, Matrix<float, 3, 1> *MeanB, Matrix<float, 6, 6> *SigA, Matrix<float, 6, 6> *SigB)
{
    Matrix<float, 4, 4, Eigen::RowMajor> X_candidate[8];
    Matrix<float, 4, 4, Eigen::RowMajor> Y_candidate[8];

    // Reshape A and B for matching the input sizes of eigen functions
    int a1 = A->rows(), a2 = A->cols(), a3 = 8;
    Map<Matrix<float, 4, 32, Eigen::RowMajor>> A_mapped(A->data(), a1, a2 * a3);
    Map<Matrix<float, 4, 32, Eigen::RowMajor>> B_mapped(B->data(), a1, a2 * a3);

    // Calculate the mean and covariance of A and B
    // non-mex function version
    // distibutionPropsMex(A_mex, &MeanA, &SigA);
    // distibutionPropsMex(B_mex, &MeanB, &SigB);
    *MeanA = A_mapped.topRows(3).rowwise().mean();
    *MeanB = B_mapped.topRows(3).rowwise().mean();
    Matrix<float, 6, 6> A_diff, B_diff;
    for (int i = 0; i < 8; ++i)
    {
        Matrix<float, 4, 4, Eigen::RowMajor> A_i = A_mapped.block(0, i * 4, 4, 4);
        Matrix<float, 4, 4, Eigen::RowMajor> B_i = B_mapped.block(0, i * 4, 4, 4);
        A_diff = A_i.topRows(3).colwise() - (*MeanA);
        B_diff = B_i.topRows(3).colwise() - (*MeanB);
        (*SigA).topLeftCorner(3, 3) += A_diff * A_diff.transpose();
        (*SigB).topLeftCorner(3, 3) += B_diff * B_diff.transpose();
        (*SigA).topRightCorner(3, 3) += A_diff * B_diff.transpose();
    }
    (*SigA) /= 8.0f;
    (*SigB) /= 8.0f;

    // update SigA and SigB if nstd_A and nstd_B are known
    if (opt)
    {
        (*SigA).topLeftCorner(3, 3) -= nstd_A * Matrix<float, 3, 3>::Identity();
        (*SigB).topLeftCorner(3, 3) -= nstd_B * Matrix<float, 3, 3>::Identity();
    }

    // Calculate Rx_solved
    EigenSolver<Matrix<float, 3, 3>> solver_A((*SigA).topLeftCorner(3, 3));
    Matrix<float
