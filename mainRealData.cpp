/*
DESCRIPTION:

The code is a MATLAB script that tests for iterative refinement on real data for
 solving the AXB=YCZ calibration problem. The problem is to find three unknown
 matrices X, Y and Z given three sets of known matrices A, B and C that satisfy
 the equation AXB=YCZ. The script uses three different methods to solve the
 problem: Prob 1, Iterative Refinement and Wang Method. Prob 1 is a probabilistic
 method that uses random sampling and voting to find an approximate solution.
 Iterative Refinement is a method that improves the solution by minimizing a
 cost function using linear least squares. Wang Method is a traditional method
 that uses singular value decomposition and matrix inversion to find an exact solution.

The script loads some real data from files that contain the matrices A, B and
 C for different poses of a robot head, hand and tag. The script also generates
 some initial guesses for X, Y and Z based on identity matrices or approximate
 measurements from kinematics data. The script then scrambles some of the data
 with different rates to simulate noise or outliers. The script then calls each
 of the three methods with the original and scrambled data and calculates the
 error metric for each method. The error metric is defined as the Frobenius norm
 of the difference between AXB and YCZ for all poses. The script then plots the
 error versus scramble rate for each method and compares their performance.

The script also defines some supporting functions that convert cell arrays to
 3D matrices, calculate the inverse of a 4x4 matrix in SE(3), calculate the adjoint
 matrix of a 4x4 matrix in SE(3), skew-symmetrize a 3x1 vector, calculate mean and
 covariance of varying data, scramble data with random permutations, construct M and
 b matrices for linear least squares problems, and calculate exponential maps of 4x4
 matrices in SE(3).

*/

/*
 * This file is equivalent to main_NAO_data_analysis.m from MATLAB
 *
 * At the moment - this file focuses on getting real data from the robot and solving for
 * X, Y, Z using Prob1 method for getting initial guesses and using that to solve using
 * Prob3 - iterative refinement method
 */

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "fKine.h"
#include "metric.h"
#include "scrambleData.h"
#include "axbyczProb1.h"
#include "axbyczProb3.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

// Convert data to 3d matrices
void convertCell2Mat(const std::vector<Eigen::MatrixXd> &headTf,
                     const std::vector<Eigen::MatrixXd> &handTf,
                     const std::vector<Eigen::MatrixXd> &tagTf,
                     std::vector<Eigen::MatrixXd> &A,
                     std::vector<Eigen::MatrixXd> &B,
                     std::vector<Eigen::MatrixXd> &C) {
    for (int i = 0; i < headTf.size(); ++i) {
        A[i] = headTf[i];
        B[i] = tagTf[i];
        C[i] = handTf[i];
    }
}

int main()
{
    // Load data
    // load('real_data/transform_ABC_unified.mat'); // load the experiment data
    // load('real_data/transform_ABC_unified_fixA.mat');
    // load('real_data/transform_ABC_unified_fixC.mat');

    // Generate Data
    // Initial guesses:
    // 1 for identity; 2(or 3) for approximate measurement from kinematics
    // data of the robot; 3 for results from Prob 1.

    int init_guess = 1;
    Eigen::Matrix4d X_init, Y_init, Z_init;

    if (init_guess == 1) {
        X_init = Eigen::Matrix4d::Identity();
        Y_init = Eigen::Matrix4d::Identity();
        Z_init = Eigen::Matrix4d::Identity();
    } else {
        Eigen::VectorXd qz1(6);
        qz1 << M_PI/6, M_PI/3, M_PI/4, M_PI/4, -M_PI/4, 0;
        Eigen::VectorXd qz2(6);
        qz2 << M_PI/3, M_PI/4, M_PI/3, -M_PI/4, M_PI/4, 0;
        Eigen::VectorXd qz3(6);
        qz3 << M_PI/4, M_PI/3, M_PI/3, M_PI/6, -M_PI/4, 0;
        X_init = fKine(qz1);
        Y_init = fKine(qz2);
        Z_init = fKine(qz3);
    }

    std::vector<std::vector<Eigen::MatrixXd>> headTf1;
    std::vector<std::vector<Eigen::MatrixXd>> handTf1;
    std::vector<std::vector<Eigen::MatrixXd>> tagTf1;
    std::vector<std::vector<Eigen::MatrixXd>> headTf2;
    std::vector<std::vector<Eigen::MatrixXd>> handTf2;
    std::vector<std::vector<Eigen::MatrixXd>> tagTf2;
    bool isRandPerm = true;

    // Choice of scramble rate
    std::vector<int> r = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    for (int rk = 0; rk < r.size(); ++rk) {
        // Convert cells to matrices
        std::vector<Eigen::MatrixXd> A;
        std::vector<Eigen::MatrixXd> B;
        std::vector<Eigen::MatrixXd> C;
        std::vector<Eigen::MatrixXd> Bp;
        for (int i = 0; i < headTf1.size(); ++i) {
            // Inputs for Iterative Refinement
            Eigen::MatrixXd A1, B1, C1;
            convertCell2Mat(headTf1[i], handTf1[i], tagTf1[i], A1, B1, C1);
            Eigen::MatrixXd A2, B2, C2;
            convertCell2Mat(headTf2[i], handTf2[i], tagTf2[i], A2, B2, C2);
            if (isRandPerm) {
                Eigen::MatrixXd Bp1 = scrambleData(B1, r[rk]);
                Eigen::MatrixXd Bp2 = scrambleData(B2, r[rk]);
            }
        }

        // Inputs for Prob 1
        std::vector<Eigen::MatrixXd> AA1(headTf1.size()), BB1(headTf1.size()), CC1(headTf1.size());
        std::vector<Eigen::MatrixXd> AA2(headTf2.size()), BB2(headTf2.size()), CC2(headTf2.size());
        std::vector<Eigen::Matrix4d> X_cal1, Y_cal1, Z_cal1;
        Eigen::MatrixXd BBp1, BBp2;

        convertCell2Mat(headTf1, handTf1, tagTf1, AA1, BB1, CC1);
        convertCell2Mat(headTf2, handTf2, tagTf2, AA2, BB2, CC2);

        if (isRandPerm) {
            BBp1 = scrambleData(BB1[0], r[rk]);
            BBp2 = scrambleData(BB2[0], r[rk]);
            Bp = scrambleData(B[0], r[rk]);
        }

        // Prob 1
        std::cout << "Probabilistic Method 1..." << std::endl;
        axbyczProb1(AA1[0], BBp1, CC1[0],
                    AA2[0], BBp2, CC2[0], 0, 0, 0,
                    X_cal1, Y_cal1, Z_cal1);

        // Initial guess for iterative refinement as the results from prob 1
        if (init_guess == 3) {
            X_init = X_cal1;
            Y_init = Y_cal1;
            Z_init = Z_cal1;
        }

        // Iterative Refinement
        std::cout << "Iterative Refinement..." << std::endl;
        axbyczProb3(A1, Bp1, C1,
                    A2, Bp2, C2,
                    X_init, Y_init, Z_init,
                    X_cal2, Y_cal2, Z_cal2,
                    num2[rk]);

        // Verification
        // Prob 1
        err1[rk] = metric(A1, B1, C1, X_cal1, Y_cal1, Z_cal1) +
                   metric(A2, B2, C2, X_cal1, Y_cal1, Z_cal1);

        // Iterative refinement
        err2[rk] = metric(A1, B1, C1, X_cal2, Y_cal2, Z_cal2) +
                   metric(A2, B2, C2, X_cal2, Y_cal2, Z_cal2);
    }


    // Plot error v.s. scramble rate
    plt::figure();
    int fontSize = 20;
    int lineW = 1;

    plt::plot(r,err1,"o-r",{{"linewidth",lineW}});
    plt::plot(r,err2,"d-g",{{"linewidth",lineW}});
    plt::plot(r,err3,"*-b",{{"linewidth",lineW}});

    auto lgd = plt::legend({"Prob 1","Iterative","Wang"});
    lgd.attr("fontsize") = fontSize;
    plt::xlabel("Scramble Rate (%)", {{"fontsize",fontSize}});
    plt::ylabel("Error",{ {"fontsize",fontSize}});

    plt::title("Real Data",{ {"fontsize",fontSize}});
}