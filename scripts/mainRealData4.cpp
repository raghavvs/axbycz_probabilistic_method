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
#include <ctime>
#include <vector>
#include <Eigen/Dense>
#include "fKine.h"
#include "metric.h"
#include "scrambleData.h"
#include "axbyczProb1.h"
#include "axbyczProb3.h"
#include "loadMatrices.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main()
{
    std::vector<Eigen::Matrix4d> A1, B1, C1, A2, B2, C2;

    std::string A1_files = {"data/20230418_abb_charuco_10x14/r1_tf.txt"};
    std::string B1_files = {"data/20230418_abb_charuco_10x14/c2b_tf.txt"};
    std::string C1_files = {"data/20230418_abb_charuco_10x14/r2_tf.txt"};
    std::string A2_files = {"data/20230418_abb_charuco_10x14/r1_tf.txt"};
    std::string B2_files = {"data/20230418_abb_charuco_10x14/c2b_tf.txt"};
    std::string C2_files = {"data/20230418_abb_charuco_10x14/r2_tf.txt"};

    loadMatrices(A1_files, A1);
    loadMatrices(B1_files, B1);
    loadMatrices(C1_files, C1);
    loadMatrices(A2_files, A2);
    loadMatrices(B2_files, B2);
    loadMatrices(C2_files, C2);

    // Resize the vectors to keep only the first 10 matrices
    int new_size = 20;
    A1.resize(new_size);
    B1.resize(new_size);
    C1.resize(new_size);
    A2.resize(new_size);
    B2.resize(new_size);
    C2.resize(new_size);

    int init_guess = 1;
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

    bool isRandPerm = true;

    std::vector<int> num_matrices_list = {2, 5, 10, 20, 50};
    std::vector<double> err1_list(num_matrices_list.size()), err3_list(num_matrices_list.size());
    // Choice of scramble rate
    std::vector<int> r = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    std::vector<double> err1(11), err3(11);
    for (int k = 0; k < num_matrices_list.size(); ++k) {
        // Resize the vectors based on the number of input matrices
        int num_matrices = num_matrices_list[k];
        A1.resize(num_matrices);
        B1.resize(num_matrices);
        C1.resize(num_matrices);
        A2.resize(num_matrices);
        B2.resize(num_matrices);
        C2.resize(num_matrices);

        for (int rk = 0; rk < r.size(); ++rk) {
            // Inputs for Prob 1
            std::vector<Eigen::Matrix4d> Bp1(A1.size());
            std::vector<Eigen::Matrix4d> Bp2(A2.size());
            std::vector<Eigen::Matrix4d> BBp1(A1.size());
            std::vector<Eigen::Matrix4d> BBp2(A2.size());

            for (int i = 0; i < A1.size(); ++i) {
                // Inputs for Iterative Refinement
                if (isRandPerm) {
                    Bp1[i] = scrambleData(B1, i, r[rk])[i];
                    Bp2[i] = scrambleData(B2, i, r[rk])[i];
                    BBp1[i] = scrambleData(B1, i, r[rk])[i];
                    BBp2[i] = scrambleData(B2, i, r[rk])[i];
                } else {
                    Bp1[i] = B1[i];
                    Bp2[i] = B2[i];
                    BBp1[i] = B1[i];
                    BBp2[i] = B2[i];
                }

            }

            // Prob 1
            //std::cout << "Probabilistic Method 1..." << std::endl;
            axbyczProb1(A1, BBp1, C1,
                        A2, BBp2, C2,
                        1, 0.001, 0.001,
                        X_cal1, Y_cal1, Z_cal1);

            // Save the errors in the respective lists
            err1_list[k] = err1[r.back()];
            //err3_list[k] = err3[r.back()];
        }

        // Plot error vs number of input matrices
        plt::figure();
        plt::plot(num_matrices_list, err1_list, "o-r");
        //plt::plot(num_matrices_list, err3_list, "s-b");
        plt::xlabel("Number of Datasets");
        plt::ylabel("Error");
        plt::title("Error vs Number of Datasets");
        plt::grid(true);
        plt::save("results/Error_vs_Number_of_Input_Matrices.png");
        plt::show();
    }
}