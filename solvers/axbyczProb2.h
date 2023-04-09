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

void axbyczProb2(const std::vector<Eigen::Matrix4d>& A1,
                 const std::vector<Eigen::Matrix4d>& B1,
                 const std::vector<Eigen::Matrix4d>& C1,
                 const std::vector<Eigen::Matrix4d>& A2,
                 const std::vector<Eigen::Matrix4d>& B2,
                 const std::vector<Eigen::Matrix4d>& C2,
                 const std::vector<Eigen::Matrix4d>& A3,
                 const std::vector<Eigen::Matrix4d>& B3,
                 const std::vector<Eigen::Matrix4d>& C3,
                 Eigen::Matrix4d& X_final,
                 Eigen::Matrix4d& Y_final,
                 Eigen::Matrix4d& Z_final) {

    // Get the first element from each input vector.
    Eigen::Matrix4d A1_fixed = A1[0];
    Eigen::Matrix4d C2_fixed = C2[0];
    Eigen::Matrix4d B3_fixed = B3[0];

    // Solve for Z
    std::vector<Eigen::Matrix4d> X, Z;
    Eigen::Matrix4d MeanA, MeanB, MeanC;
    Eigen::Matrix<double, 6, 6> SigA, SigB, SigC;

    std::vector<Eigen::Matrix4d> Z_g, X_dummy, Y_dummy;
    Eigen::Matrix4d MeanC1, MeanB1, MeanA2, MeanB2;
    Eigen::Matrix<double, 6, 6> SigC1, SigB1, SigA2, SigB2;

    batchSolveXY(C1, B1, false,0,0,Z_g,Y_dummy,
                 MeanC1,MeanB1,SigC1,SigB1);

    // Keep the candidates of Z that are SE3
    for (const auto& Z_candidate : Z_g) {
        if (Z_candidate.determinant() > 0) {
            Z.push_back(Z_candidate);
        }
    }

    // Solve for X
    std::vector<Eigen::Matrix4d> X_g;
    std::vector<Eigen::Matrix4d> A2_inv(A2.size()), B2_inv(B2.size());

    for (size_t i = 0; i < A2.size(); ++i) {
        A2_inv[i] = A2[i].inverse();
        B2_inv[i] = B2[i].inverse();
    }

    batchSolveXY(A2, B2_inv, false, 0, 0, X_g, Y_dummy,
                 MeanA2, MeanB2, SigB, SigC);

    batchSolveXY(A2_inv, B2, false, 0, 0, X_g, Y_dummy,
                 MeanA2, MeanB2, SigB, SigC);

    // Keep the candidates of X that are SE3
    for (const auto& X_candidate : X_g) {
        if (X_candidate.determinant() > 0) {
            X.push_back(X_candidate);
        }
    }

    // Solve for Y
    std::vector<Eigen::Matrix4d> Y_g_inv, Y;
    Eigen::Matrix4d MeanA3, MeanC3;
    std::vector<Eigen::Matrix4d> A3_inv(A3.size()), C3_inv(C3.size());

    for (size_t i = 0; i < A3.size(); ++i) {
        A3_inv[i] = A3[i].inverse();
        C3_inv[i] = C3[i].inverse();
    }

    batchSolveXY(C3_inv, A3_inv, false, 0, 0, Y_g_inv, Y_dummy,
                 MeanC, MeanA, SigC, SigA);

    batchSolveXY(C3, A3, false, 0, 0, Y_g_inv, Y_dummy,
                 MeanC3, MeanA3, SigC, SigA);

    // Keep the candidates of Y that are SE3
    for (const auto& Y_candidate_inv : Y_g_inv) {
        if (Y_candidate_inv.determinant() > 0) {
            Y.push_back(Y_candidate_inv.inverse());
        }
    }

    // Find out the optimal (X, Y, Z) that minimizes cost
    size_t s_X = X.size();
    size_t s_Y = Y.size();
    size_t s_Z = Z.size();

    Eigen::MatrixXd cost(s_X, s_Y * s_Z);
    double weight = 1.8; // weight on the translational error of the cost function
    double min_cost = std::numeric_limits<double>::max();
    int min_i = 0, min_j = 0, min_p = 0;

    for (size_t i = 0; i < s_X; ++i) {
        for (size_t j = 0; j < s_Y; ++j) {
            for (size_t p = 0; p < s_Z; ++p) {
                Eigen::Matrix4d left1 = A1_fixed * X[i] * MeanB1;
                Eigen::Matrix4d right1 = Y[j] * MeanC1 * Z[p];
                double diff1 = rotError(left1, right1) + weight * tranError(left1, right1);

                Eigen::Matrix4d left2 = MeanA2 * X[i] * MeanB2;
                Eigen::Matrix4d right2 = Y[j] * C2_fixed * Z[p];
                double diff2 = rotError(left2, right2) + weight * tranError(left2, right2);

                Eigen::Matrix4d left3 = MeanA3 * X[i] * B3_fixed;
                Eigen::Matrix4d right3 = Y[j] * MeanC3 * Z[p];
                double diff3 = rotError(left3, right3) + weight * tranError(left3, right3);

                double current_cost = std::abs(diff1) + std::abs(diff2) + std::abs(diff3);
                cost(i, j * s_Z + p) = current_cost;
                if (current_cost < min_cost) {
                    min_cost = current_cost;
                    min_i = i;
                    min_j = j;
                    min_p = p;
                }
            }
        }
    }

    //// Recover the X, Y, Z that minimizes cost
    X_final = X[min_i];
    Z_final = Z[min_j];
    Y_final = Y[min_p];
}

#endif