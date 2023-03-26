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
#include "permFixABC.h"
#include "generateSetsOfABC.h"
#include "axbyczProb1.h"
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

    for (int s = 1; s <= num; ++s) {
        int counter = 0;
        int optPDF = 1;
        Eigen::MatrixXd A1, B1, C1, A2, B2, C2, A3, B3, C3;
        std::tie(A1, B1, C1, A2, B2, C2, A3, B3, C3) = generateSetsOfABC(Num, optPDF, gmean, k*Cov, XActual, YActual, ZActual);
        for (double r : rate) {
            counter++;
            Eigen::MatrixXd A1_perm, B1_perm, C1_perm;
            std::tie(A1_perm,B1_perm,C1_perm) = permFixABC(A1.block(0,0,A1.rows(),A1.cols()), B1,C1,r);
            Eigen::MatrixXd C2_perm,A2_perm,B2_perm;
            std::tie(C2_perm,A2_perm,B2_perm) = permFixABC(C2.block(0,0,C2.rows(),C2.cols()), A2,B2,r);
            Eigen::MatrixXd B3_perm,C3_perm,A3_perm;
            std::tie(B3_perm,C3_perm,A3_perm) = permFixABC(B3.block(0,0,B3.rows(),B3.cols()), C3,A3,r);
            Eigen::MatrixXd A_perm(A1.rows(),A1.cols(),A1.size()+A2.size()+A3.size());
            A_perm << A1_perm,A2_perm,A3_perm;
            Eigen::MatrixXd B_perm(B1.rows(),B1.cols(),B1.size()+B2.size()+B3.size());
            B_perm << B1_perm,B2_perm,B3_perm;
            Eigen::MatrixXd C_perm(C1.rows(),C1.cols(),C1.size()+C2.size()+C3.size());
            C_perm << C1_perm,C2_perm,C3_perm;

            // different algorithms
            axbyczProb1(A1, B1_perm, C1_perm, A2_perm, B2_perm, C2, 0, 0, 0, X_f1, Y_f1, Z_f1);
            axbyczProb2(A1, B1_perm, C1_perm, A2_perm, B2_perm, C2, A3_perm, B3, C3, X_f2, Y_f2, Z_f2);

            // Error analysis
            // Prob1 Error
            Err_1.block(counter,0,1,3) = getErrorAXBYCZ(X_f1,Y_f1,Z_f1,XActual,YActual,ZActual).transpose();
            // Prob2 Error
            Err_2.block(counter,0,1,3) = getErrorAXBYCZ(X_f2,Y_f2,Z_f2,XActual,YActual,ZActual).transpose();

        }
    }

    // Compute the averaged errors
    Eigen::MatrixXd Err1_Avg = Err_1.colwise().sum() / num;
    Eigen::MatrixXd Err2_Avg = Err_2.colwise().sum() / num;
    Eigen::MatrixXd ErrW_Avg = Err_W.colwise().sum() / num;
    Eigen::MatrixXd ErrNP_Avg = Err_NP.colwise().sum() / num;
    Eigen::MatrixXd ErrDK_Avg = Err_DK.colwise().sum() / num;

    // Plots
    plt::figure();
    std::vector<std::string> y_lb = {"$E_{\\bf R_{X}}$", "$E_{\\bf R_{Y}}$", "$E_{\\bf R_{Z}}$",
                                     "$E_{\\bf t_{X}}$", "$E_{\\bf t_{Y}}$", "$E_{\\bf t_{Z}}$"};

    for (int i = 0; i < 6; ++i) {
        // Subplot 1
        plt::subplot(2,3,i+1);
        plt::plot(rate, Err1_Avg.col(i), "d-r");
        plt::plot(rate, Err2_Avg.col(i), "*-g");

        plt::legend({"Prob1", "Prob2"});
        plt::ylabel(y_lb[i]);
    }
    plt::show();

    return 0;
}

