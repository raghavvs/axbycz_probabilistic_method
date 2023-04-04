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
#include "axbyczProb3.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
    int init_guess = 2;
    Eigen::Matrix4d X_init, Y_init, Z_init;

    if (init_guess == 1) {
        X_init.setIdentity();
        Y_init.setIdentity();
        Z_init.setIdentity();
    }

    bool isRandPerm = true;
    std::vector<int> r = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    int num_trials = 20;

    Eigen::VectorXd err_prob, err_iter;
    double err1, err2;
    //Eigen::MatrixXd err1_mat, err2_mat;
    Eigen::MatrixXd err1_mat = Eigen::MatrixXd::Zero(r.size(), num_trials);
    Eigen::MatrixXd err2_mat = Eigen::MatrixXd::Zero(r.size(), num_trials);

    for (int n = 0; n < num_trials; ++n) {
        Eigen::Matrix4d X_true, Y_true, Z_true;
        initializeXYZ(1, X_true, Y_true, Z_true);

        for (int rk = 0; rk < r.size(); ++rk) {
            int length = 100;
            Eigen::VectorXd Mean = Eigen::VectorXd::Zero(6);
            Eigen::MatrixXd Cov = 0.1 * Eigen::MatrixXd::Identity(6, 6);

            std::vector<Eigen::Matrix4d> A, B, C, Bp;

            std::vector<Eigen::Matrix4d> A1, B1, C1, A2, B2, C2;

            for (int i = 0; i < 5; ++i) {
                std::tie(A1, B1, C1) = generateABC(length, 1, 1, Mean, Cov,
                                                                X_true, Y_true, Z_true);
                std::tie(A2, B2, C2) = generateABC(length, 3, 1, Mean, Cov,
                                                                X_true, Y_true, Z_true);

                A.insert(A.end(), A1.begin(), A1.end());
                A.insert(A.end(), A2.begin(), A2.end());

                B.insert(B.end(), B1.begin(), B1.end());
                B.insert(B.end(), B2.begin(), B2.end());

                C.insert(C.end(), C1.begin(), C1.end());
                C.insert(C.end(), C2.begin(), C2.end());

                if (isRandPerm) {
                    Bp.emplace_back(scrambleData(B1[i], r[rk]));
                    Bp.emplace_back(scrambleData(B2[i], r[rk]));
                }
            }

            // Inputs for Prob 1
            Eigen::Matrix4d AA1 = A1[0];
            Eigen::Matrix4d BB1 = B1[0];
            Eigen::Matrix4d CC1 = C1[0];
            Eigen::Matrix4d AA2 = A2[0];
            Eigen::Matrix4d BB2 = B2[0];
            Eigen::Matrix4d CC2 = C2[0];
            Eigen::MatrixXd BBp1, BBp2;

            if (isRandPerm) {
                BBp1 = scrambleData(BB1, r[rk]);
                BBp2 = scrambleData(BB2, r[rk]);
                Bp.emplace_back(scrambleData(BB1, r[rk]));
            }

            // Prob 1
            Eigen::Matrix4d X_cal1, Y_cal1, Z_cal1;

            axbyczProb1(AA1, BBp1, CC1, AA2, BBp2, CC2,
                        true, 0.001, 0.001,
                        X_cal1, Y_cal1, Z_cal1);

            if (init_guess == 2) {
                X_init = X_cal1;
                Y_init = Y_cal1;
                Z_init = Z_cal1;
            }


            Eigen::Matrix4d X_cal2, Y_cal2, Z_cal2;
            int num2;

            axbyczProb3(A1, B1, C1, A2, B2, C2,
                        X_init, Y_init, Z_init,
                        X_cal2, Y_cal2, Z_cal2, num2);

            // Verification
            // Prob 1
            err_prob = getErrorAXBYCZ(X_cal1, Y_cal1, Z_cal1,
                                      X_true, Y_true, Z_true);
            err1 = metric(A1, B1, C1, X_cal1, Y_cal1, Z_cal1)
                            + metric(A2, B2, C2, X_cal1, Y_cal1, Z_cal1);
            err1_mat.coeffRef(rk, n) = err1;

            // Iterative refinement
            err_iter = getErrorAXBYCZ(X_cal2, Y_cal2, Z_cal2,
                                                      X_true, Y_true, Z_true);
            err2 = metric(A1, B1, C1, X_cal2, Y_cal2, Z_cal2)
                            + metric(A2, B2, C2, X_cal2, Y_cal2, Z_cal2);
            err2_mat.coeffRef(rk, n) = err2;
        }
    }

    // Compute the averaged errors
    Eigen::MatrixXd err_prob_avg = err_prob.colwise().sum() / num_trials;
    Eigen::MatrixXd err_iter_avg = err_iter.colwise().sum() / num_trials;

    Eigen::VectorXd err1_avg = err1_mat.rowwise().sum() / num_trials;
    Eigen::VectorXd err2_avg = err2_mat.rowwise().sum() / num_trials;

    std::cout << "works till here? - YES" << std::endl;

    std::cout << "err1_avg: " << err1_avg << std::endl;
    std::cout << "err2_avg: " << err2_avg << std::endl;

    std::cout << "works till here?" << std::endl;

     ///// PLOTS

     // Plot error v.s. scramble rate
    // Errors with ground truth
    plt::figure();
    int fontSize = 20;
    int lineW = 1;
    std::vector<std::string> y_lb = {"$R_{X}$", "$R_{Y}$", "$R_{Z}$",
                                     "${\\bf t_{X}}$", "${\\bf t_{Y}}$", "${\\bf t_{Z}}$"};

    auto x_data = Eigen::Map<const Eigen::VectorXd>(reinterpret_cast<const double *>(r.data()), r.size());

    for (int i = 0; i < 6; ++i) { //

        plt::subplot(2, 3, i+1);

        // Plot data
        auto y_data = err_prob_avg.col(i).cast<double>();
        std::vector<double> x_data_vec(x_data.data(), x_data.data() + x_data.size());
        std::vector<double> y_data_vec(y_data.data(), y_data.data() + y_data.size());
        std::cout << "Size of x_data_vec: " << x_data_vec.size() << ", Size of y_data_vec: " << y_data_vec.size() << std::endl;
        plt::plot(x_data_vec, y_data_vec, "o-r");

        // Set axis labels
        plt::xlabel("Scramble Rate / %");
        plt::ylabel(y_lb[i]);

        // Set legend
        plt::legend();
    }

    // Errors between two sides of calibration equations
    plt::figure();
    std::vector<double> r_vec(r.begin(), r.end());
    std::vector<double> err1_avg_vec(err1_avg.data(), err1_avg.data() + err1_avg.size());
    std::vector<double> err2_avg_vec(err2_avg.data(), err2_avg.data() + err2_avg.size());
    plt::plot(r_vec, err1_avg_vec, "o-r");
    plt::plot(r_vec, err2_avg_vec, "d-g");

    plt::legend();
    plt::xlabel("Scramble Rate (%)");
    plt::ylabel("Error");
    plt::show();
}
