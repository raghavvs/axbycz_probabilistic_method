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
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "fKine.h"
#include "metric.h"
#include "scrambleData.h"
#include "axbyczProb1.h"
#include "axbyczProb2.h"
#include "axbyczProb3.h"
#include "loadMatrices.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main()
{
    std::vector<Eigen::Matrix4d> A1, B1, C1, A2, B2, C2;

    std::vector<std::string> A1_files = {"data/006.txt"};
    std::vector<std::string> B1_files = {"data/006_board.txt"};
    std::vector<std::string> C1_files = {"data/006_m.txt"};
    std::vector<std::string> A2_files = {"data/007.txt"};
    std::vector<std::string> B2_files = {"data/007_board.txt"};
    std::vector<std::string> C2_files = {"data/007_m.txt"};

    loadMatrices(A1_files, A1);
    loadMatrices(B1_files, B1);
    loadMatrices(C1_files, C1);
    loadMatrices(A2_files, A2);
    loadMatrices(B2_files, B2);
    loadMatrices(C2_files, C2);

    // Generate Data
    // Initial guesses:
    // 1 for identity; 2(or 3) for approximate measurement from kinematics
    // data of the robot; 3 for results from Prob 1.

    int init_guess = 2;
    Eigen::Matrix4d X_init, Y_init, Z_init;
    Eigen::Matrix4d X_cal1, Y_cal1, Z_cal1, X_cal3, Y_cal3, Z_cal3;

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

    bool isRandPerm = false;

    // Choice of scramble rate
    std::vector<int> r = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    std::vector<double> err1(11), err2(11), err3(11);
    for (int rk = 0; rk < r.size(); ++rk) {
        std::vector<Eigen::MatrixXd> A;
        std::vector<Eigen::MatrixXd> B;
        std::vector<Eigen::MatrixXd> C;
        for (int i = 0; i < A1.size(); ++i) {
            // Inputs for Iterative Refinement
            if (isRandPerm) {
                Eigen::MatrixXd Bp1 = scrambleData(B1[i], r[rk]);
                Eigen::MatrixXd Bp2 = scrambleData(B2[i], r[rk]);
            }
        }

        // Inputs for Prob 1
        std::vector<Eigen::Matrix4d> AA1(A1.size()), BB1(A1.size()), CC1(A1.size());
        std::vector<Eigen::Matrix4d> AA2(A2.size()), BB2(A2.size()), CC2(A2.size());
        Eigen::Matrix4d BBp1, BBp2, Bp;

        if (isRandPerm) {
            BBp1 = scrambleData(BB1[0], r[rk]);
            BBp2 = scrambleData(BB2[0], r[rk]);
            Bp = scrambleData(B[0], r[rk]);
        }

        // Prob 1
        std::cout << "Probabilistic Method 1..." << std::endl;
        /*axbyczProb1(AA1[0], BBp1, CC1[0],
                    AA2[0], BBp2, CC2[0],
                    0, 0, 0,
                    X_cal1, Y_cal1, Z_cal1);*/

        axbyczProb1(A1[0], B1[0], C1[0],
                    A2[0], B2[0], C2[0],
                    0, 0, 0,
                    X_cal1, Y_cal1, Z_cal1);


        // Prob 2
        std::cout << "Probabilistic Method 2..." << std::endl;
        axbyczProb2(A1[0], B1[0], C1[0],
                    A2[0], B2[0], C2[0],
                    A3[0], B3[0], C3[0],
                    X_cal2, Y_cal2, Z_cal2);



        // Initial guess for iterative refinement as the results from prob 1
        if (init_guess == 3) {
            X_init = X_cal1;
            Y_init = Y_cal1;
            Z_init = Z_cal1;
        }

        // Iterative Refinement
        std::cout << "Iterative Refinement..." << std::endl;
        int num2[rk];
        std::vector<Eigen::Matrix4d> Bp1, Bp2;
        axbyczProb3(A1, B1, C1,
                    A2, B2, C2,
                    X_init, Y_init, Z_init,
                    X_cal3, Y_cal3, Z_cal3,
                    num2[rk]);

        /*axbyczProb3(A1, Bp1, C1,
                    A2, Bp2, C2,
                    X_init, Y_init, Z_init,
                    X_cal2, Y_cal2, Z_cal2,
                    num2[rk]);*/

        // Verification
        // Prob 1
        err1[rk] = metric(A1, B1, C1, X_cal1, Y_cal1, Z_cal1) +
                            metric(A2, B2, C2, X_cal1, Y_cal1, Z_cal1);
        // Prob 2
        err2[rk] = metric(A1, B1, C1, X_cal2, Y_cal2, Z_cal2) +
                   metric(A2, B2, C2, X_cal2, Y_cal2, Z_cal2) +
                   metric(A3, B3, C3, X_cal2, Y_cal2, Z_cal2);

        // Iterative refinement
        err3[rk] = metric(A1, B1, C1, X_cal3, Y_cal3, Z_cal3) +
                   metric(A2, B2, C2, X_cal3, Y_cal3,Z_cal3);
    }

    std::ofstream outFile("results/XYZ.txt", std::ios_base::app);

    outFile << "Probability Method 1" << std::endl;
    outFile << "X_cal1: " << std::endl << X_cal1 << std::endl;
    outFile << "Y_cal1: " << std::endl << Y_cal1 << std::endl;
    outFile << "Z_cal1: " << std::endl << Z_cal1 << std::endl;
    outFile << "Probability Method 2" << std::endl;
    outFile << "X_cal2: " << std::endl << X_cal2 << std::endl;
    outFile << "Y_cal2: " << std::endl << Y_cal2 << std::endl;
    outFile << "Z_cal2: " << std::endl << Z_cal2 << std::endl;
    outFile << "Probability Method 3 - Iterative Refinement" << std::endl;
    outFile << "X_cal3: " << std::endl << X_cal3 << std::endl;
    outFile << "Y_cal3: " << std::endl << Y_cal3 << std::endl;
    outFile << "Z_cal3: " << std::endl << Z_cal3 << std::endl;

    outFile.close();

    // Plot error vs scramble rate
    plt::figure();

    plt::plot(r, err1, "o-r");
    plt::plot(r, err2, "d-g");
    plt::plot(r, err3, "s-b");

    // Choose appropriate x and y coordinates for the labels
    double label_x = r.back() * 1.05;
    double label_y1 = err1.back();
    double label_y2 = err2.back();
    double label_y3 = err3.back();

    plt::text(label_x, label_y1, "Method 1");
    plt::text(label_x, label_y2, "Method 2");
    plt::text(label_x, label_y3, "Method 3");

    plt::xlabel("Scramble Rate (%)");
    plt::ylabel("Error");

    plt::title("Real Data");

    plt::grid(true);

    plt::save("results/Error_vs_Scramble_Rate_10.png");

    plt::show();
}
