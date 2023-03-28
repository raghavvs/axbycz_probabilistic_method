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

Input:
    A1, B1, C1, A2, B2, C2: Matrices - dim 4x4
    opt: bool
    nstd1, nst2: standard deviation
Output:
    X_final, Y_final, Z_final: Matrices - dim 4x4
*/

#ifndef AXBYCZPROB1_H
#define AXBYCZPROB1_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include "batchSolveXY.h"
#include "rotError.h"
#include "tranError.h"

void axbyczProb1(const Eigen::Matrix4d &A1,
                 const Eigen::Matrix4d &B1,
                 const Eigen::Matrix4d &C1,
                 const Eigen::Matrix4d &A2,
                 const Eigen::Matrix4d &B2,
                 const Eigen::Matrix4d &C2,
                 bool opt,
                 double nstd1,
                 double nstd2,
                 std::vector<Eigen::Matrix4d> &X_final,
                 std::vector<Eigen::Matrix4d> &Y_final,
                 std::vector<Eigen::Matrix4d> &Z_final) {

    //   A1 is constant with B1 and C1 free
    //   C2 is constant with A2 and B2 free

    int len = 8;
    std::vector<Eigen::Matrix4d> Z_g(len);
    Eigen::MatrixXd MeanA, MeanB, MeanC, SigA, SigB, SigC;

    std::vector<Eigen::Matrix4d> A1_vec(len), B1_vec(len), C1_vec(len),
            A2_vec(len), B2_vec(len), C2_vec(len);
    for (int i = 0; i < len; ++i) {
        A1_vec[i] = A1.block(0, 0, 4, 4);
        B1_vec[i] = B1.block(0, 0, 4, 4);
        C1_vec[i] = C1.block(0, 0, 4, 4);
        A2_vec[i] = A2.block(0, 0, 4, 4);
        B2_vec[i] = B2.block(0, 0, 4, 4);
        C2_vec[i] = C2.block(0, 0, 4, 4);
    }

    //// ------ using probability methods ------
    // calculate Z_g : all guesses of Z
    //// ------ Solve for Z -------- //
    // A1 fixed, B1 and C1 free

    Eigen::MatrixXd MeanB1, MeanC1, MeanA2;

    batchSolveXY(C1_vec, B1_vec, len, opt,nstd1,nstd2,Z_g,Y_final,
                 MeanC1,MeanB1,SigC,SigB);

    // Keep the candidates of Z that are SE3
    // Normally there will be four Z \in SE3
    int Z_index = 0;
    for (int i = 0; i < Z_g.size(); ++i) {
        if (Z_g[i].determinant() > 0) {
            Z_final.push_back(Z_g[i]);
            ++Z_index;
        }
    }

    int s_Z = Z_final.size();

    //// ------ Solve for X -------- //
    // C2 fixed, A2 and B2 free

    // ------ Calculate B2^-1 -------
    int Num = A2.size();
    std::vector<Eigen::Matrix4d> A2_inv(Num), B2_inv(Num);
    for (int i = 0; i < Num; ++i) {
        A2_inv[i] = A2_vec[i].inverse();
        B2_inv[i] = B2_vec[i].inverse();
    }

    // ------ using probability methods ------
    // calculate X_g : all guesses of X
    std::vector<Eigen::Matrix4d> X_g(len);
    batchSolveXY(A2_vec, B2_inv, len, opt, nstd1, nstd2, X_g, Y_final,
                 MeanA2, MeanB, SigC, SigB);

    // Calculate MeanB for computing Y later
    // Note: can be further simplified by using only the distribution function
    Eigen::MatrixXd MeanB2;
    batchSolveXY(A2_inv, B2_vec, len, opt, nstd1, nstd2, X_g, Y_final,
                 MeanA, MeanB2, SigC, SigB);

    // Keep the candidates of X that are SE3å
    // Normally there will be four X \in SE3
    int X_index = 0;
    std::vector<Eigen::Matrix4d> X;
    for (int i = 0; i < X_g.size(); ++i) {
        if (X_g[i].determinant() > 0) {
            X.push_back(X_g[i]);
            ++X_index;
        }
    }

    int s_X = X.size();

    //// ------ Solve for Y -------- //
    // Compute Y using the mean equations
    size_t dim = 2 * s_X * s_Z;

    std::vector<Eigen::Matrix4d> Y(dim, Eigen::Matrix4d::Zero());
    for (int i = 0; i < s_X; ++i) {
        for (int j = 0; j < s_Z; ++j) {
            Y[i * s_Z + j] = (A1 * X[i] * MeanB1 * Z_final[j].inverse()) * MeanC1.inverse();
            Y[i * s_Z + j + s_X * s_Z] = (MeanA2 * X[i] * MeanB2) * Z_final[j].inverse() * C2.inverse();
        }
    }

    int s_Y = Y.size();

    //// Find out the optimal (X, Y, Z) that minimizes cost

    Eigen::MatrixXd cost = Eigen::MatrixXd::Zero(s_X, s_Y * s_Z);
    double weight = 1.5; // weight on the translational error of the cost function

    for (int i = 0; i < s_X; ++i) {
        for (int j = 0; j < s_Z; ++j) {
            for (int m = 0; m < s_Y; ++m) {
                Eigen::MatrixXd left1 = A1 * X[i] * MeanB1;
                Eigen::MatrixXd right1 = Y[m] * MeanC1 * Z_final[j];

                double diff1 =
                        rotError(left1, right1) + weight * tranError(left1, right1);

                Eigen::MatrixXd left2 = MeanA2 * X[i] * MeanB2;
                Eigen::MatrixXd right2 = Y[m] * C2 * Z_final[j];
                double diff2 =
                        rotError(left2, right2) + weight * tranError(left2, right2);

                // different error metrics can be picked and this (diff1 +
                // diff2) is the best one so far. However, it can still be
                // unstable sometimes and miss the optimal solutions
                cost(i, j * s_Y + m) = std::norm(diff1) + std::norm(diff2);
            }
        }
    }

    //// recover the X,Y,Z that minimizes cost

    Eigen::MatrixXd::Index minRow, minCol;
    cost.minCoeff(&minRow, &minCol);
    int I1 = minRow * cost.cols() + minCol;

    int I_row = I1 / cost.cols();
    int I_col = I1 % cost.cols();

    Eigen::Matrix4d X_final_ = X[I_row]; // final X

    int index_Z;
    if (I_col % s_Y > 0) {
        index_Z = floor(I_col/s_Y) + 1;
    } else {
        index_Z = floor(I_col/s_Y);
    }

    Eigen::Matrix4d Z_final_ = Z_final[index_Z]; // final Z

    int index_Y;
    if (I_col % s_Y > 0) {
        index_Y = I_col % s_Y;
    } else {
        index_Y = s_Y;
    }

    Eigen::Matrix4d Y_final_ = Y[index_Y]; // final Y

    std::cout << "X_final_: " << std::endl << X_final_ << std::endl;
    std::cout << "Y_final_: " << std::endl << Y_final_ << std::endl;
    std::cout << "Z_final_: " << std::endl << Z_final_ << std::endl;
}

#endif