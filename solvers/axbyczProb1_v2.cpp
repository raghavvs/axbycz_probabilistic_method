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
#include <vector>
#include <eigen3/Eigen/Dense>
#include "batchSolveXY.h"

void axbyczProb1(const Eigen::Matrix4d& A1,
                 const Eigen::Matrix4d& B1,
                 const Eigen::Matrix4d& C1,
                 const Eigen::Matrix4d& A2,
                 const Eigen::Matrix4d& B2,
                 const Eigen::Matrix4d& C2,
                 bool opt,
                 double nstd1,
                 double nstd2,
                 int len,
                 std::vector<Eigen::MatrixXd>& X_final,
                 std::vector<Eigen::MatrixXd>& Y_final,
                 std::vector<Eigen::MatrixXd>& Z_final) {

    //   A1 is constant with B1 and C1 free
    //   C2 is constant with A2 and B2 free

    //// ------ Solve for Z -------- //
    // A1 fixed, B1 and C1 free

    //// ------ using probability methods ------

    std::vector<Eigen::MatrixXd> MeanA(2), MeanB(2), MeanC(2), SigA(2), SigB(2), SigC(2);

    std::vector<Eigen::MatrixXd> Z_g;

    batchSolveXY(C1, B1, len, opt, nstd1, nstd2, X_final, Z_final,
                 MeanC[0], MeanB[0], SigC[0], SigB[0]);

    std::vector<Eigen::MatrixXd> Z;

    for (const auto &z: Z_g) {
        if (z.determinant() > 0) {
            Z.push_back(z);
        }
    }

    int s_Z = static_cast<int>(Z.size());

    std::cout << "works till here?" << std::endl;

    //// ------ Solve for X -------- //
    // C2 fixed, A2 and B2 free

    //// ------ Calculate B2^-1 -------

    size_t dim1 = 4;
    size_t dim2 = 4;
    int Num = static_cast<int>(A2.size());

    Eigen::MatrixXd A2_inv, B2_inv;

    for (int i = 0; i < Num; ++i) {
        Eigen::Matrix4d A2_mat = A2.block(0, i*dim2, dim1, dim2);
        Eigen::Matrix4d B2_mat = B2.block(0, i*dim2, dim1, dim2);
        A2_inv(i) = A2_mat.inverse();
        B2_inv(i) = B2_mat.inverse();
    }

    //// ------ using probability methods ------
    // calculate X_g : all guesses of X
    std::vector<Eigen::MatrixXd> X_g, MeanA2, MeanB2;

    batchSolveXY(A2, B2, len, opt, nstd1, nstd2, X_g, Y_final,
                 MeanA2[0], MeanB2[0], SigA[0], SigB[0]);

    // Calculate MeanB2 for computing Y later
    // Note: can be further simplified by using only the distribution function

    batchSolveXY(A2_inv, B2, len, opt, nstd1, nstd2, X_final, Y_final,
                 MeanA2[0], MeanB2[0], SigA[0], SigB[0]);

    // Keep the candidates of X that are SE3
    // Normally, there will be four X \in SE3
    int X_index = 1;

    std::vector<Eigen::MatrixXd> X;

    for (const auto &x : X_g) {
        if (x.determinant() > 0) {
            X.push_back(x);
            ++X_index;
        }
    }

}