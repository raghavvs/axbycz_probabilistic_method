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

void axbyczProb1(const std::vector<Eigen::Matrix4d>& A1,
                 const std::vector<Eigen::Matrix4d>& B1,
                 const std::vector<Eigen::Matrix4d>& C1,
                 const std::vector<Eigen::Matrix4d>& A2,
                 const std::vector<Eigen::Matrix4d>& B2,
                 const std::vector<Eigen::Matrix4d>& C2,
                 bool opt,
                 double nstd1,
                 double nstd2,
                 Eigen::Matrix4d& X_final,
                 Eigen::Matrix4d& Y_final,
                 Eigen::Matrix4d& Z_final) {

    //std::cout << "axbyczProb1:" << std::endl;

    ////   A1 is constant with B1 and C1 free
    Eigen::Matrix4d A1_fixed = A1[0];

    ////   C2 is constant with A2 and B2 free
    Eigen::Matrix4d C2_fixed = C2[0];

    // Solve for Z
    std::vector<Eigen::Matrix4d> Z_g, X_dummy, Y_dummy;
    Eigen::Matrix4d MeanC1, MeanB1, MeanA2, MeanB2;
    Eigen::Matrix<double, 6, 6> SigC1, SigB1, SigA2, SigB2;

    batchSolveXY(C1, B1, opt, nstd1, nstd2, Z_g, Y_dummy,
                 MeanC1, MeanB1, SigC1, SigB1);

    std::vector<Eigen::Matrix4d> Z;
    for (const auto& z : Z_g) {
        if (z.determinant() > 0) {
            Z.push_back(z);
        }
    }

    size_t s_Z = Z.size();

    // Calculate B2_inv
    int Num = A2.size();
    std::vector<Eigen::Matrix4d> A2_inv(Num), B2_inv(Num);
    for (int i = 0; i < Num; ++i) {
        A2_inv[i] = A2[i].inverse();
        B2_inv[i] = B2[i].inverse();
    }

    // Solve for X
    std::vector<Eigen::Matrix4d> X_g;
    batchSolveXY(A2, B2_inv, opt, nstd1, nstd2, X_g, Y_dummy,
                 MeanA2, MeanB2, SigA2, SigB2);

    std::vector<Eigen::Matrix4d> X;
    for (const auto& x : X_g) {
        if (x.determinant() > 0) {
            X.push_back(x);
        }
    }

    size_t s_X = X.size();

    // Calculate MeanB2 for computing Y later
    batchSolveXY(A2_inv, B2, opt, nstd1, nstd2, X_dummy, Y_dummy,
                 MeanA2, MeanB2, SigA2, SigB2);

    // Compute Y
    std::vector<Eigen::Matrix4d> Y(2 * s_X * s_Z);
    for (size_t i = 0; i < s_X; ++i) {
        for (size_t j = 0; j < s_Z; ++j){
            Eigen::Matrix4d left = A1_fixed * X[i] * MeanB1;
            Eigen::Matrix4d right = Z[j].eval().inverse() * MeanC1.eval().inverse();
            Y[(i * s_Z) + j] = left * right;

            left = MeanA2 * X[i] * MeanB2;
            right = C2_fixed * Z[j].eval().inverse();
            Y[(i * s_Z) + j + s_X * s_Z] = left * right;
        }
    }

    size_t s_Y = Y.size();

    // Find the optimal (X, Y, Z) that minimizes cost
    Eigen::MatrixXd cost(s_X, s_Y * s_Z);
    double weight = 1.5;
    double min_cost = std::numeric_limits<double>::max();
    int min_i = 0, min_j = 0, min_m = 0;

    for (size_t i = 0; i < s_X; ++i) {
        for (size_t j = 0; j < s_Z; ++j) {
            for (size_t m = 0; m < s_Y; ++m) {
                Eigen::Matrix4d left1 = A1_fixed * X[i] * MeanB1;
                Eigen::Matrix4d right1 = Y[m] * MeanC1 * Z[j];

                double diff1 = rotError(left1, right1) + weight *
                                                         tranError(left1, right1);

                Eigen::Matrix4d left2 = MeanA2 * X[i] * MeanB2;
                Eigen::Matrix4d right2 = Y[m] * C2_fixed * Z[j];

                double diff2 = rotError(left2, right2) + weight *
                                                         tranError(left2, right2);

                double current_cost = std::abs(diff1) + std::abs(diff2);
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

    //// Recover the X, Y, Z that minimize cost
    X_final = X[min_i];
    Z_final = Z[min_j];
    Y_final = Y[min_m];

}

#endif