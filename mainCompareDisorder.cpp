/*
DESCRIPTION:

This code initializes several variables and arrays such as counter, gmean, Cov, k, num, Num, and rate.
 It also calls a function named initializeXYZ with an input argument of 2 and assigns its outputs to
 the variables XActual, YActual, and ZActual.

The code then initializes several error containers as 3D arrays of zeros with dimensions based on the
 length of the variable rate, 6, and the variable num.

The main part of the code is a nested for loop that iterates over the values in the variable num and
 rate. Within this loop, it calls a function named generateSetsOfABC with several input arguments
 including Num, optPDF (which is set to 1), gmean, k*Cov, XActual, YActual.

 Within the inner for loop, the code calls several functions with various input arguments to
 permute and concatenate data streams. It then calls several different algorithms to solve for
 X, Y, and Z values using the permuted and concatenated data. These algorithms include axbyczProb1,
 axbyczProb2, Wang_AXBYCZ, Yan_AXBYCZ_PN, and Yan_AXBYCZ_DK.

The code then performs error analysis by calling a function named getErrorAXBYCZ with several
 input arguments including the calculated values of X, Y, and Z from the different algorithms
 as well as the actual values of XActual, YActual, and ZActual. The resulting errors are
 stored in the previously initialized error containers.

Overall, it appears that this code is performing simulations to compare the performance of
 different algorithms for solving a problem involving data streams A, B, and C.

*/

/*
 * This file is equivalent to main_comparison_disorder.m from MATLAB
 * It compares the error for varying scramble rate
 */

#include <iostream>
#include <Eigen/Dense>
#include "initializeXYZ.h"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

int main() {
    int counter = 0;
    Eigen::VectorXd gmean(6);
    gmean << 0, 0, 0, 0, 0, 0;
    Eigen::MatrixXd Cov = Eigen::MatrixXd::Identity(6, 6);
    double k = 0.1;
    int num = 1; // number of simulations
    int Num = 80; // number of data
    Eigen::VectorXi rate(6);
    rate << 0, 20, 40, 60, 80, 100;

    // generate random X, Y and Z
    int opt_XYZ = 2;
    Eigen::Matrix4d XActual, YActual, ZActual;
    initializeXYZ(opt_XYZ, XActual, YActual, ZActual);

    // Error containers initialization
    std::vector<std::vector<std::vector<double>>> Err_1(rate.size(), std::vector<std::vector<double>>(6, std::vector<double>(num)));
    std::vector<std::vector<std::vector<double>>> Err_2(rate.size(), std::vector<std::vector<double>>(6, std::vector<double>(num)));
    std::vector<std::vector<std::vector<double>>> Err_W(rate.size(), std::vector<std::vector<double>>(6, std::vector<double>(num)));
    std::vector<std::vector<std::vector<double>>> Err_NP(rate.size(), std::vector<std::vector<double>>(6, std::vector<double>(num)));
    std::vector<std::vector<std::vector<double>>> Err_DK(rate.size(), std::vector<std::vector<double>>(6, std::vector<double>(num)));

    for (int s = 1; s <= num; s++) {
        counter = 0;
        int optPDF = 1;
        auto [A1, B1, C1, A2, B2, C2, A3, B3, C3] =
                generateSetsOfABC(Num, optPDF, gmean, k * Cov, XActual, YActual, ZActual);

        for (auto r: rate) {
            counter++;
            auto [A1_perm, B1_perm, C1_perm] = permFixABC(
                    A1.slice(std::array<long long int, 3>{0LL, 0LL, (long long int) (s - 1)}), B1, C1, r);
            auto [C2_perm, A2_perm, B2_perm] = permFixABC(
                    C2.slice(std::array<long long int, 3>{0LL, 0LL, (long long int) (s - 1)}), A2, B2, r);
            auto [B3_perm, C3_perm, A3_perm] = permFixABC(
                    B3.slice(std::array<long long int, 3>{0LL, 0LL, (long long int) (s - 1)}), C3, A3, r);

            // concatenate data for traditional simultaneous axbycz solvers
            std::vector<Eigen::MatrixXd> A_vec{A1_perm, A2_perm, A3_perm};
            std::vector<Eigen::MatrixXd> B_vec{B1_perm, B2_perm, B3_perm};
            std::vector<Eigen::MatrixXd> C_vec{C1_perm, C2_perm, C3_perm};

            Eigen::Tensor<double, 3> A_perm(A1.rows(), A1.cols(), A_vec.size());
            Eigen::Tensor<double, 3> B_perm(B1.rows(), B1.cols(), B_vec.size());
            Eigen::Tensor<double, 3> C_perm(C1.rows(), C1.cols(), C_vec.size());

            for (int i = 0; i < A_vec.size(); i++) {
                A_perm.chip(i, 2) = A_vec[i];
                B_perm.chip(i, 2) = B_vec[i];
                C_perm.chip(i, 2) = C_vec[i];
            }

            // different algorithms
            auto [X_f1, Y_f1, Z_f1] = axbyczProb1(A1, B1_perm, C1_perm, A2_perm, B2_perm, C2, 0, 0, 0);
            auto [X_f2, Y_f2, Z_f2] = axbyczProb2(A1, B1_perm, C1_perm, A2_perm, B2_perm, C2, A3_perm, B3, C3);
            auto [X_wang, Y_wang, Z_wang] = Wang_AXBYCZ(A.perm, B.perm, C.perm, XActual, YActual, ZActual);
            auto [X_NP, Y_NP, Z_NP] = Yan_AXBYCZ_PN(A.perm, B.perm, C.perm, XActual, YActual, ZActual);
            auto [X_DK, Y_DK, Z_DK] = Yan_AXBYCZ_DK(A1, B1_perm, C1_perm, A2_perm, B2_perm, C2);

            // err analysis
            for (int counter = 0; counter < rate.size(); counter++) {
                // Prob1 Error
                Err_1.chip(counter, 0).chip(s - 1, 2) = getErrorAXBYCZ(X_f1, Y_f1, Z_f1, XActual, YActual, ZActual);
                // Prob2 Error
                Err_2.chip(counter, 0).chip(s - 1, 2) = getErrorAXBYCZ(X_f2, Y_f2, Z_f2, XActual, YActual, ZActual);
                // Wang Error
                Err_W.chip(counter, 0).chip(s - 1, 2) = getErrorAXBYCZ(X_wang, Y_wang, Z_wang, XActual, YActual,
                                                                       ZActual);
                // NP Error
                Err_NP.chip(counter, 0).chip(s - 1, 2) = getErrorAXBYCZ(X_NP, Y_NP, Z_NP, XActual, YActual, ZActual);
                // DK Error
                Err_DK.chip(counter, 0).chip(s - 1, 2) = getErrorAXBYCZ(X_DK, Y_DK, Z_DK, XActual, YActual, ZActual);
            }

            // Compute the averaged errors
            Eigen::MatrixXd Err1_Avg = Err_1.sum(Eigen::array<int, 1>{2}) / num;
            Eigen::MatrixXd Err2_Avg = Err_2.sum(Eigen::array<int, 1>{2}) / num;
            Eigen::MatrixXd ErrW_Avg = Err_W.sum(Eigen::array<int, 1>{2}) / num;
            Eigen::MatrixXd ErrNP_Avg = Err_NP.sum(Eigen::array<int, 1>{2}) / num;
            Eigen::MatrixXd ErrDK_Avg = Err_DK.sum(Eigen::array<int, 1>{2}) / num;
        }
    }

    // Plots
    plt::figure();
    std::vector<std::string> y_lb{"$E_{\\bf R_{X}}$", "$E_{\\bf R_{Y}}$", "$E_{\\bf R_{Z}}$",
                                  "$E_{\\bf t_{X}}$", "$E_{\\bf t_{Y}}$", "$E_{\\bf t_{Z}}$"};

    for (int i = 0; i < 6; i++) {
        // Subplot 1
        plt::subplot(2, 3, i + 1);
        plt::plot(rate, Err1_Avg.col(i), "d-r");
        plt::plot(rate, Err2_Avg.col(i), "*-g");
        plt::plot(rate, ErrW_Avg.col(i), "o-b");
        plt::plot(rate, ErrNP_Avg.col(i), "+-m");
        plt::plot(rate, ErrDK_Avg.col(i), "x-c");

        auto len1 = plt::legend({"Prob1", "Prob2", "Wang", "PN", "DK"});
        len1.attr("fontsize") = 14;
        len1.attr("interpreter") = "latex";
        plt::ylabel(y_lb[i], {{"fontsize",    18},
                              {"interpreter", "latex"}});
    }
    plt::show();

    return 0;
}

