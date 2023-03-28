/*
DESCRIPTION:

The script generates data by setting initial guesses for X, Y, and Z to
 either identity or results from Prob 1 based on the value of init_guess.
 The script then sets the choice of scramble rate and the number of trials.
 For each trial, the script displays the current trial number and generates
 true values for simulations using the InitializeXYZ function. For each value
 of r, the script sets the length, mean, and covariance matrix for data
 generation. The script then generates A, B, and C using different data
 distributions and noise levels by calling the generateABC function. If
 isRandPerm is true, the script scrambles the data in B1 and B2 using the
 scrambleData function. The inputs for Prob 1 are then set and the script
 displays that it is running Probabilistic Method 1.

*/

/*
 * This file is equivalent to main_NAO_simulation.m from MATLAB
 */

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "initializeXYZ.h"
#include "generateABC.h"
#include "scrambleData.h"
#include "getErrorAXBYCZ.h"
#include "metric.h"
#include "axbyczProb1.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
    // Generate Data
    // Initial guesses:
    // 1 for identity; 2 for results from Prob 1.
    int init_guess = 2;
    Eigen::Matrix4d X_init, Y_init, Z_init;
    // Initial guess as Identity
    if (init_guess == 1) {
        X_init = Eigen::Matrix4d::Identity();
        Y_init = Eigen::Matrix4d::Identity();
        Z_init = Eigen::Matrix4d::Identity();
    }

    // Choice of scramble rate
    bool isRandPerm = true;
    std::vector<int> r = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    int num_trials = 20;

    for (int n = 1; n <= num_trials; ++n) {
        std::cout << "Num of trials: " << n << std::endl;
        // True values for simulations
        Eigen::Matrix4d X_true, Y_true, Z_true;
        std::tie(X_true, Y_true, Z_true) = initializeXYZ(1);

        for (int rk = 0; rk < r.size(); ++rk) {
            int length = 100;
            Eigen::VectorXd Mean(6);
            Mean << 0, 0, 0, 0, 0, 0;
            Eigen::MatrixXd Cov = 0.1 * Eigen::MatrixXd::Identity(6, 6);

            std::vector<Eigen::MatrixXd> A, B, C, Bp;
            std::vector<std::vector<Eigen::MatrixXd>> A1(5), B1(5), C1(5), A2(5), B2(5), C2(5);

            // Generate A, B, C using different data distribution and noise level
            for (int i = 0; i < 5; ++i) {
                std::tie(A1[i], B1[i], C1[i]) = generateABC(length, 1, 1, Mean, Cov, X_true, Y_true, Z_true);
                std::tie(A2[i], B2[i], C2[i]) = generateABC(length, 3, 1, Mean, Cov, X_true, Y_true, Z_true);

                Eigen::MatrixXd AA(A1[i].rows(), A1[i].cols(), A1[i].size() + A2[i].size());
                AA << A1[i], A2[i];
                Eigen::MatrixXd BB(B1[i].rows(), B1[i].cols(), B1[i].size() + B2[i].size());
                BB << B1[i], B2[i];
                Eigen::MatrixXd CC(C1[i].rows(), C1[i].cols(), C1[i].size() + C2[i].size());
                CC << C1[i], C2[i];

                A.push_back(AA);
                B.push_back(BB);
                C.push_back(CC);

                if (isRandPerm) {
                    Eigen::MatrixXd Bp1 = scrambleData(B1[i], r[rk]);
                    Eigen::MatrixXd Bp2 = scrambleData(B2[i], r[rk]);
                }
            }

            // Inputs for Prob 1
            Eigen::MatrixXd AA1 = A1[0];
            Eigen::MatrixXd BB1 = B1[0];
            Eigen::MatrixXd CC1 = C1[0];
            Eigen::MatrixXd AA2 = A2[0];
            Eigen::MatrixXd BB2 = B2[0];
            Eigen::MatrixXd CC2 = C2[0];
            Eigen::MatrixXd BBp1, BBp2, Bp;
            if (isRandPerm) {
                BBp1 = scrambleData(BB1, r[rk]);
                BBp2 = scrambleData(BB2, r[rk]);
                Bp = scrambleData(B, r[rk]);
            }

            // Prob 1
            std::cout << "Probabilistic Method 1..." << std::endl;
            Eigen::Matrix4d X_cal1, Y_cal1, Z_cal1;
            std::tie(X_cal1, Y_cal1, Z_cal1) = axbyczProb1(AA1.block(0, 0, AA1.rows(), AA1.cols()), BBp1,
                                                           CC1.block(0, 0, CC1.rows(), CC1.cols()), AA2, BBp2,
                                                           CC2.block(0, 0, CC2.rows(), CC2.cols()), 0, 0, 0);

            // Initial guess for iterative refinement as the results from prob 1
            if (init_guess == 2) {
                X_init = X_cal1;
                Y_init = Y_cal1;
                Z_init = Z_cal1;
            }

            // Iterative Refinement
            std::cout << "Iterative Refinement..." << std::endl;
            Eigen::Matrix4d X_cal2, Y_cal2, Z_cal2;
            int num2;
            std::tie(X_cal2, Y_cal2, Z_cal2, num2) = axbyczProb3(A1, Bp1, C1, A2, Bp2, C2, X_init, Y_init, Z_init);

            // Verification
            // Prob 1
            Eigen::Vector3d err_prob = getErrorAXBYCZ(X_cal1, Y_cal1, Z_cal1, X_true, Y_true, Z_true);
            double err1 = metric(A1, B1, C1, X_cal1, Y_cal1, Z_cal1) + metric(A2, B2, C2, X_cal1, Y_cal1, Z_cal1);

            // Iterative refinement
            Eigen::Vector3d err_iter = getErrorAXBYCZ(X_cal2, Y_cal2, Z_cal2, X_true, Y_true, Z_true);
            double err2 = metric(A1, B1, C1, X_cal2, Y_cal2, Z_cal2) + metric(A2, B2, C2, X_cal2, Y_cal2, Z_cal2);
        }
    }

    // Compute the averaged errors
    Eigen::MatrixXd err_prob_avg = err_prob.colwise().sum() / num_trials;
    Eigen::MatrixXd err_iter_avg = err_iter.colwise().sum() / num_trials;
    Eigen::MatrixXd err_wang_avg = err_wang.colwise().sum() / num_trials;

    Eigen::VectorXd err1_avg = err1.rowwise().sum() / num_trials;
    Eigen::VectorXd err2_avg = err2.rowwise().sum() / num_trials;
    Eigen::VectorXd err3_avg = err3.rowwise().sum() / num_trials;

    // Plot error v.s. scramble rate
    // Errors with ground truth
    plt::figure();
    int fontSize = 20;
    int lineW = 1;
    std::vector<std::string> y_lb = {"$R_{X}$", "$R_{Y}$", "$R_{Z}$",
                                     "${\\bf t_{X}}$", "${\\bf t_{Y}}$", "${\\bf t_{Z}}$"};

    for (int i = 0; i < 6; ++i) {
        // Subplot 1
        plt::subplot(2, 3, i + 1);
        plt::plot(r, err_prob_avg.col(i), "o-r", {{"linewidth", lineW}});
        plt::plot(r, err_iter_avg.col(i), "d-g", {{"linewidth", lineW}});
        plt::plot(r, err_wang_avg.col(i), "*-b", {{"linewidth", lineW}});

        plt::legend({"Prob1", "Iterative", "Wang"});
        plt::ylabel(y_lb[i]);
        plt::xlabel("Scramble Rate / %");
    }
    // Errors between two sides of calibration equations
    plt::figure();
    plt::plot(r, err1_avg, "o-r", {{"linewidth", lineW}});
    plt::plot(r, err2_avg, "d-g", {{"linewidth", lineW}});
    plt::plot(r, err3_avg, "*-b", {{"linewidth", lineW}});

    plt::legend({"Prob 1", "Iterative", "Wang"});
    plt::xlabel("Scramble Rate (%)");
    plt::ylabel("Error");
    plt::show();
}
