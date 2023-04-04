/*
DESCRIPTION:

This function implements the Prob2 method in the paper
Prerequisites on the input
  A1 is constant with B1 and C1 free
  C2 is constant with A1 adn B1 free
  B3 is constant with A3 and C3 free

Input:
    A1, B1, C1, A2, B2, C2: Matrices - dim 4x4
    opt: bool
    nstd1, nst2: standard deviation
Output:
    X_final, Y_final, Z_final: Matrices - dim 4x4
*/

#ifndef AXBYCZPROB2_H
#define AXBYCZPROB2_H

#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include "batchSolveXY.h"
#include "rotError.h"
#include "tranError.h"

void axbyczProb2(const Eigen::Matrix4d &A1,
                 const Eigen::Matrix4d &B1,
                 const Eigen::Matrix4d &C1,
                 const Eigen::Matrix4d &A2,
                 const Eigen::Matrix4d &B2,
                 const Eigen::Matrix4d &C2,
                 const Eigen::Matrix4d &A3,
                 const Eigen::Matrix4d &B3,
                 const Eigen::Matrix4d &C3,
                 Eigen::Matrix4d &X_final,
                 Eigen::Matrix4d &Y_final,
                 Eigen::Matrix4d &Z_final) {

    // A1 fixed, B1 and C1 free
    int len = 8;
    double nstd1 = 0;
    double nstd2 = 0;
    bool opt = false;

    std::vector<Eigen::Matrix4d> Z_g(len), X(len), Y_temp(len), Z(len);
    Eigen::MatrixXd MeanA, MeanB, MeanC, SigA, SigB, SigC;

    std::vector<Eigen::Matrix4d> A1_vec(len), B1_vec(len), C1_vec(len),
            A2_vec(len), B2_vec(len), C2_vec(len),
            A3_vec(len), B3_vec(len), C3_vec(len);
    for (int i = 0; i < len; ++i) {
        A1_vec[i] = A1.block(0, 0, 4, 4);
        B1_vec[i] = B1.block(0, 0, 4, 4);
        C1_vec[i] = C1.block(0, 0, 4, 4);
        A2_vec[i] = A2.block(0, 0, 4, 4);
        B2_vec[i] = B2.block(0, 0, 4, 4);
        C2_vec[i] = C2.block(0, 0, 4, 4);
        A3_vec[i] = A3.block(0, 0, 4, 4);
        B3_vec[i] = B3.block(0, 0, 4, 4);
        C3_vec[i] = C3.block(0, 0, 4, 4);
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

    //// ------ Solve for Y -------- //
    // B3 fixed, A3 and C3 free

    // ------ Calculate B2^-1 -------
    std::vector<Eigen::Matrix4d> A3_inv(len);
    std::vector<Eigen::Matrix4d> C3_inv(len);
    for (int i = 0; i < len; i++) {
        A3_inv[i] = A3_vec[i].inverse();
        C3_inv[i] = C3_vec[i].inverse();
    }

    // ------ using probability methods ------
    // calculate X_g : all guesses of X
    std::vector<Eigen::Matrix4d> Y_g_inv;
    batchSolveXY(C3_inv, A3_inv, len, opt, nstd1, nstd2, Y_g_inv, Y_temp,
                 MeanC, MeanA, SigC, SigA);

    // Calculate MeanA2 and MeanC2 for the cost function later
    Eigen::MatrixXd MeanC3, MeanA3;
    batchSolveXY(C3_vec, A3_vec, len, opt, nstd1, nstd2, Y_g_inv, Y_temp,
                 MeanC3, MeanA3, SigC, SigA);

    // Keep the candidates of Y that are SE3
    int Y_index = 1;
    std::vector<Eigen::Matrix4d> Y(4, Eigen::Matrix4d::Zero());
    for (int i = 0; i < Y_g_inv.size(); i++) {
        if (Y_g_inv[i].determinant() > 0) {
            Y[Y_index] = Y_g_inv[i].inverse();
            Y_index++;
        }
    }

    //// Find out the optimal (X, Y, Z) that minimizes cost

    int s_Z = Z_final.size();
    int s_X = X.size();
    int s_Y = Y.size();


    Eigen::MatrixXd cost = Eigen::MatrixXd::Zero(s_X, s_Y * s_Z);
    double weight = 1.8; // weight on the translational error of the cost function
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