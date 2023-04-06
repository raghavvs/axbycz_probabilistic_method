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

In the case of two robotic arms:
     A - robot 1's base to end effector transformation (forward kinematics)
     B - camera to calibration target transformation
     C - robot 2's base to end effector transformation (forward kinematics)
     X - end effector of robot 1 to camera transformation
     Y - robot 1's base to robot 2's base transformation
     Z - end effector of robot 2 to calibration target transformation
*/


#ifndef AXBYCZPROB1_H
#define AXBYCZPROB1_H

#include <iostream>
#include <vector>
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
                 Eigen::Matrix4d &X_final,
                 Eigen::Matrix4d &Y_final,
                 Eigen::Matrix4d &Z_final) {

    //   A1 is constant with B1 and C1 free
    //   C2 is constant with A2 and B2 free

    int len = 8;
    std::vector<Eigen::Matrix4d> Z_g(len), X(len), Y_temp(len), Z(len);
    Eigen::MatrixXd MeanA, MeanB, MeanC, SigA, SigB, SigC;

    std::vector<Eigen::Matrix4d> A1_vec(len), B1_vec(len), C1_vec(len),
                                    A2_vec(len), B2_vec(len), C2_vec(len);
    for (int i = 0; i < len; ++i) {
        A1_vec[i] = A1;
        B1_vec[i] = B1;
        C1_vec[i] = C1;
        A2_vec[i] = A2;
        B2_vec[i] = B2;
        C2_vec[i] = C2;
    }

    //// ------ using probability methods ------
    // calculate Z_g : all guesses of Z
    //// ------ Solve for Z -------- //
    // A1 fixed, B1 and C1 free

    Eigen::MatrixXd MeanB1, MeanC1, MeanA2;

    batchSolveXY(C1_vec, B1_vec, len, opt,nstd1,nstd2,Z_g,Y_temp,
                 MeanC1,MeanB1,SigC,SigB);

    // Keep the candidates of Z that are SE3
    // Normally there will be four Z \in SE3
    int Z_index = 0;
    for (int i = 0; i < Z_g.size(); ++i) {
        if (Z_g[i].determinant() > 0) {
            Z.push_back(Z_g[i]);
            ++Z_index;
        }
    }

    int s_Z = Z.size();

    //// ------ Solve for X -------- //
    // C2 fixed, A2 and B2 free

    // ------ Calculate B2^-1 -------
    int Num = A2_vec.size();
    std::vector<Eigen::Matrix4d> A2_inv(Num), B2_inv(Num);
    for (int i = 0; i < Num; ++i) {
        A2_inv[i] = A2_vec[i].inverse();
        B2_inv[i] = B2_vec[i].inverse();
    }

    // ------ using probability methods ------
    // calculate X_g : all guesses of X
    std::vector<Eigen::Matrix4d> X_g(len);
    batchSolveXY(A2_vec, B2_inv, len, opt, nstd1, nstd2, X_g, Y_temp,
                 MeanA2, MeanB, SigC, SigB);

    // Calculate MeanB for computing Y later
    // Note: can be further simplified by using only the distribution function
    Eigen::MatrixXd MeanB2;
    batchSolveXY(A2_inv, B2_vec, len, opt, nstd1, nstd2, X_g, Y_temp,
                 MeanA, MeanB2, SigC, SigB);

    // Keep the candidates of X that are SE3å
    // Normally there will be four X \in SE3
    int X_index = 0;
    for (int i = 0; i < X_g.size(); ++i) {
        if (X_g[i].determinant() > 0) {
            X.push_back(X_g[i]);
            ++X_index;
        }
    }

    int s_X = X.size();

    //// ------ Solve for Y -------- //
    // Compute Y using the mean equations
    std::vector<Eigen::Matrix4d> Y(2*s_X*s_Z);
    for (int i = 0; i < s_X; ++i) {
        for (int j = 0; j < s_Z; ++j) {
            Y[i * s_Z + j] = (A1 * X[i] * MeanB1 * Z[j].inverse()) * MeanC1.inverse();
            Y[i * s_Z + j + s_X * s_Z] = (MeanA2 * X[i] * MeanB2) * Z[j].inverse() * C2.inverse();
        }
    }

    int s_Y = Y.size();

    //// Find out the optimal (X, Y, Z) that minimizes cost

    Eigen::MatrixXd cost = Eigen::MatrixXd::Zero(s_X, s_Y * s_Z);
    double weight = 1.5; // weight on the translational error of the cost function
    double min_cost = std::numeric_limits<double>::max();
    int min_i = 0, min_j = 0, min_m = 0;

    for (int i = 0; i < s_X; ++i) {
        for (int j = 0; j < s_Z; ++j) {
            for (int m = 0; m < s_Y; ++m) {
                Eigen::MatrixXd left1 = A1 * X[i] * MeanB1;
                Eigen::MatrixXd right1 = Y[m] * MeanC1 * Z[j];

                double diff1 =
                        rotError(left1, right1) + weight * tranError(left1, right1);

                Eigen::MatrixXd left2 = MeanA2 * X[i] * MeanB2;
                Eigen::MatrixXd right2 = Y[m] * C2 * Z[j];
                double diff2 =
                        rotError(left2, right2) + weight * tranError(left2, right2);

                double current_cost = diff1 + diff2;
                cost(i, j * s_Y + m) = current_cost;

                if (current_cost < min_cost) {
                    min_cost = current_cost;
                    min_i = i;
                    min_j = j;
                    min_m = m;
                }
            }
        }
    }

    //// recover the X,Y,Z that minimizes cost
    X_final = X[min_i];
    Z_final = Z[min_j];
    Y_final = Y[min_m];
}

#endif