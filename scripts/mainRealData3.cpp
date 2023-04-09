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
#include <ctime>
#include <vector>
#include <Eigen/Dense>
#include "fKine.h"
#include "metric.h"
#include "scrambleData.h"
#include "axbyczProb3.h"
#include "loadMatrices.h"
#include "matplotlibcpp.h"
#include <cmath>
#include "batchSolveXY.h"
#include "rotError.h"
#include "tranError.h"

namespace plt = matplotlibcpp;

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

    ////   A1 is constant with B1 and C1 free
    Eigen::Matrix4d A1_fixed = A1[0];

    ////   C2 is constant with A2 and B2 free
    Eigen::Matrix4d C2_fixed = C2[0];

    std::cout << A1_fixed << std::endl;
    std::cout << C2_fixed << std::endl;

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

    std::cout << X.size() << std::endl;
    std::cout << Y.size() << std::endl;
    std::cout << Z.size() << std::endl;

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

    std::cout << cost.rows() << cost.cols() << std::endl;

    //// Recover the X, Y, Z that minimize cost
    X_final = X[min_i];
    Z_final = Z[min_j];
    Y_final = Y[min_m];
}

int main()
{
    std::vector<Eigen::Matrix4d> A1, B1, C1, A2, B2, C2;

    std::vector<std::string> A1_files = {"data/r1_tf1.txt"};
    std::vector<std::string> B1_files = {"data/c2b_tf1.txt"};
    std::vector<std::string> C1_files = {"data/r2_tf1.txt"};
    std::vector<std::string> A2_files = {"data/r1_tf1.txt"};
    std::vector<std::string> B2_files = {"data/c2b_tf1.txt"};
    std::vector<std::string> C2_files = {"data/r2_tf1.txt"};

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
    Eigen::Matrix4d X_cal1, Y_cal1, Z_cal1, X_cal2, Y_cal2, Z_cal2, X_cal3, Y_cal3, Z_cal3;

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
    std::vector<double> err1(11), err3(11);
    for (int rk = 0; rk < 10; ++rk) {
        std::vector<Eigen::Matrix4d> Bp1, Bp2, BBp1, BBp2;

        // Inputs for Prob 1
        Bp1.resize(A1.size());
        Bp2.resize(A2.size());
        BBp1.resize(A1.size());
        BBp2.resize(A2.size());

        for (int i = 0; i < A1.size(); ++i) {
            // Inputs for Iterative Refinement
            if (isRandPerm) {
                Bp1[i] = scrambleData(B1[i], r[rk]);
                Bp2[i] = scrambleData(B2[i], r[rk]);
                BBp1[i] = scrambleData(B1[i], r[rk]);
                BBp2[i] = scrambleData(B2[i], r[rk]);
            }
        }

        // Prob 1
        //std::cout << "Probabilistic Method 1..." << std::endl;
        /*axbyczProb1(A1, BBp1, C1,
                    A2, BBp2, C2,
                    1, 0.0001, 0.0001,
                    X_cal1, Y_cal1, Z_cal1);*/
        axbyczProb1(A1, B1, C1,
                    A2, B2, C2,
                    1, 0.0001, 0.0001,
                    X_cal1, Y_cal1, Z_cal1);

        // Initial guess for iterative refinement as the results from prob 1
        if (init_guess == 3) {
            X_init = X_cal1;
            Y_init = Y_cal1;
            Z_init = Z_cal1;
        }

        // Iterative Refinement
        //std::cout << "Iterative Refinement..." << std::endl;
        int num = 1;

        /*axbyczProb3(A1, Bp1, C1,
                    A2, Bp2, C2,
                    X_init, Y_init, Z_init,
                    X_cal3, Y_cal3, Z_cal3 ,
                    num);*/
        axbyczProb3(A1, B1, C1,
                    A2, B2, C2,
                    X_init, Y_init, Z_init,
                    X_cal3, Y_cal3, Z_cal3 ,
                    num);

        // Verification
        // Prob 1
        err1[rk] = metric(A1, B1, C1, X_cal1, Y_cal1, Z_cal1) +
                            metric(A2, B2, C2, X_cal1, Y_cal1, Z_cal1);

        // Iterative refinement
        err3[rk] = metric(A1, B1, C1, X_cal3, Y_cal3, Z_cal3) +
                   metric(A2, B2, C2, X_cal3, Y_cal3,Z_cal3);
    }

    std::ofstream outFile("results/XYZ.txt", std::ios_base::app);

    // Get current date and time
    std::time_t now = std::time(nullptr);
    char buffer[100];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    outFile << "Current date and time: " << buffer << std::endl;

    outFile << "Probability Method 1" << std::endl;
    outFile << "X_cal1: " << std::endl << X_cal1 << std::endl;
    outFile << "Y_cal1: " << std::endl << Y_cal1 << std::endl;
    outFile << "Z_cal1: " << std::endl << Z_cal1 << std::endl;

    outFile << "Probability Method 3 - Iterative Refinement" << std::endl;
    outFile << "X_cal3: " << std::endl << X_cal3 << std::endl;
    outFile << "Y_cal3: " << std::endl << Y_cal3 << std::endl;
    outFile << "Z_cal3: " << std::endl << Z_cal3 << std::endl;

    outFile.close();

    std::cout << "Error 1 [0]: " << err1[0] << std::endl;
    std::cout << "Error 1 [10]: " << err1[100] << std::endl;
    std::cout << "Error 3 [0]: " << err3[0] << std::endl;
    std::cout << "Error 3 [10]: " << err3[100] << std::endl;

   /* // Plot error vs scramble rate
    plt::figure();

    plt::plot(r, err1, "o-r");
    plt::plot(r, err3, "s-b");

    // Choose appropriate x and y coordinates for the labels
    double label_x = r.back() * 1.05;
    double label_y1 = err1.back();
    double label_y3 = err3.back();

    plt::text(label_x, label_y1, "Method 1");
    plt::text(label_x, label_y3, "Method 3");

    plt::xlabel("Scramble Rate (%)");
    plt::ylabel("Error");

    plt::title("Real Data");

    plt::grid(true);

    plt::save("results/Error_vs_Scramble_Rate_16.png");

    plt::show();*/
}