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
#include "axbyczProb2.h"
#include "matplotlibcpp.h"
#include "plotProbResults.h"

namespace plt = matplotlibcpp;

int main() {
    // Initialize Parameters
    int counter = 0;
    Eigen::VectorXd gmean(6);
    gmean << 0, 0, 0, 0, 0, 0;
    std::vector<double> coeff1{0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14}; // Scaling factor for the covariances
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(6,6); // coeff1*cov
    std::vector<double> coeff2{0.02, 0.06, 0.10, 0.14}; // Scaling factor for the covariances
    Eigen::MatrixXd cov_noise = Eigen::MatrixXd::Identity(6,6); // coeff2*cov
    int num = 10; // number of simulations
    int Num = 100; // number of data
    std::string optPlot = "lineplot"; // Plot the averaged error : 'lineplot' & ''boxplot'

    int opt_XYZ = 2; // generate random X, Y and Z
    Eigen::Matrix4d XActual, YActual, ZActual;
    initializeXYZ(opt_XYZ, XActual, YActual, ZActual);

    // Error container initialization
    if (optPlot == "boxplot") {
        std::vector<Eigen::MatrixXd> Err11(coeff1.size(), Eigen::MatrixXd::Zero(num,6));
        std::vector<Eigen::MatrixXd> Err21(coeff1.size(), Eigen::MatrixXd::Zero(num,6));
    } else if (optPlot == "lineplot") {
        std::vector<Eigen::MatrixXd> Err11(coeff1.size(), Eigen::MatrixXd::Zero(6, num));
        std::vector<Eigen::MatrixXd> Err21(coeff1.size(), Eigen::MatrixXd::Zero(6, num));
        std::vector<Eigen::MatrixXd> Err12(coeff1.size(), Eigen::MatrixXd::Zero(6, num));
        std::vector<Eigen::MatrixXd> Err22(coeff1.size(), Eigen::MatrixXd::Zero(6, num));
    }

    Eigen::MatrixXd A11, B11, C11, A21, B21, C21, A31, B31, C31,
                    A12, B12, C12, A22, B22, C22, A32, B32, C32,
                    A13, B13, C13, A23, B23, C23, A33, B33, C33;

    std::vector<Eigen::Matrix4d> X_f11, Y_f11, Z_f11,
                                X_f12, Y_f12, Z_f12,
                                X_f13, Y_f13, Z_f13,
                                X_f21, Y_f21, Z_f21,
                                X_f22, Y_f22, Z_f22,
                                X_f23, Y_f23, Z_f23;

    std::vector<Eigen::MatrixXd> Err_1, Err_2;

    for (double k : coeff1) {
        counter++;
        for (int s = 1; s <= num; s++) {
            // Generate data triples with different distributions
            int optPDF = 1;
            std::tie(A11, B11, C11, A21, B21, C21, A31, B31, C31) = generateSetsOfABC(Num, optPDF, gmean, k * Cov,
                                                                                      XActual, YActual, ZActual);

            optPDF = 2;
            std::tie(A12, B12, C12, A22, B22, C22, A32, B32, C32) = generateSetsOfABC(Num, optPDF, gmean, k * Cov,
                                                                                      XActual, YActual, ZActual);

            optPDF = 3;
            std::tie(A13, B13, C13, A23, B23, C23, A33, B33, C33) = generateSetsOfABC(Num, optPDF, gmean, k * Cov,
                                                                                      XActual, YActual, ZActual);

            // Solve for X, Y and Z
            axbyczProb1(A11, B11, C11, A21, B21, C21, 0, 0, 0, X_f11, Y_f11, Z_f11);
            axbyczProb2(A11, B11, C11, A21, B21, C21, A31, B31, C31, X_f21, Y_f21, Z_f21);

            axbyczProb1(A12, B12, C12, A22, B22, C22, A32, B32, C32, 0, 0, 0, X_f12, Y_f12, Z_f12);
            axbyczProb2(A11, B11, C11, A21, B21, C21, A31, B31, C31, X_f22, Y_f22, Z_f22);

            axbyczProb1(A13, B13, C13, A23, B23, C23, A33, B33, C33, 0, 0, 0, X_f13, Y_f13, Z_f13);
            axbyczProb2(A11, B11, C11, A21, B21, C21, A31, B31, C31, X_f23, Y_f23, Z_f23);

            // Error Analysis
            if (optPlot == "boxplot") {
                // Prob1 Error
                Err11[counter - 1].row(s - 1) = getErrorAXBYCZ(X_f1, Y_f1, Z_f1, XActual, YActual, ZActual);
            } else if (optPlot == "lineplot") {
                // Prob1 Error with Data of 1st Distribution
                Err11[counter - 1].col(s - 1) = getErrorAXBYCZ(X_f11, Y_f11, Z_f11, XActual, YActual, ZActual);
                // Prob1 Error with Data of 2nd Distribution
                Err12[counter - 1].col(s - 1) = getErrorAXBYCZ(X_f12, Y_f12, Z_f12, XActual, YActual, ZActual);
                // Prob1 Error with Data of 3rd Distribution
                Err13[counter - 1].col(s - 1) = getErrorAXBYCZ(X_f13, Y_f13, Z_f13, XActual, YActual, ZActual);
                // Prob2 Error with Data of 1st Distribution
                Err21[counter - 1].col(s - 1) = getErrorAXBYCZ(X_f21, Y_f21, Z_f21, XActual, YActual, ZActual);
                // Prob2 Error with Data of 2nd Distribution
                Err22[counter - 1].col(s - 1) = getErrorAXBYCZ(X_f22, Y_f22, Z_f22, XActual, YActual, ZActual);
                // Prob2 Error with Data of 3rd Distribution
                Err23[counter - 1].col(s - 1) = getErrorAXBYCZ(X_f23, Y_f23, Z_f23, XActual, YActual, ZActual);
            }

            // Plot X and Y
            if (false) {
                try {
                    plt::figure((counter - 1) * 3 + 1);
                    // trplot(XActual, 'color', 'r');
                    // plt::hold_on();
                    // trplot(X_f11,'color','b');
                    // plt::axis("auto");
                    // plt::hold_on();
                    plt::title("Final X");

                    plt::figure((counter - 1) * 3 + 2);
                    // trplot(YActual, 'color', 'r');
                    // plt::hold_on();
                    // trplot(Y_f11,'color','b');
                    // plt::axis("auto");
                    // plt::hold_on();
                    plt::title("Final Y");

                    plt::figure((counter - 1) * 3 + 3);
                    // trplot(ZActual, 'color', 'r');
                    // plt::hold_on();
                    // trplot(Z_f11,'color','b');
                    // plt::axis("auto");
                    // plt::hold_on();
                    plt::title("Final Z");
                } catch (...) {
                    std::cout << Y << std::endl;
                }
            }
        }
    }

    // Plot the averaged error w.r.t. covariances
    plotProbResults(Err11, Err21, coeff1, optPlot);
    plotProbResults(Err12, Err22, coeff1, optPlot);
    plotProbResults(Err13, Err23, coeff1, optPlot);
}