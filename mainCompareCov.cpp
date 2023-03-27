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
    // [XActual,YActual,ZActual] = InitializeXYZ(opt_XYZ);

    // Error container initialization
    if (optPlot == "boxplot") {
        std::vector<Eigen::MatrixXd> Err11(coeff1.size(), Eigen::MatrixXd::Zero(num,6));
        std::vector<Eigen::MatrixXd> Err21(coeff1.size(), Eigen::MatrixXd::Zero(num,6));
    } else if (optPlot == "lineplot") {
        std::vector<Eigen::MatrixXd> Err11(coeff1.size(), Eigen::MatrixXd::Zero(6,num));
        std::vector<Eigen::MatrixXd> Err21(coeff1.size(), Eigen::MatrixXd::Zero(6,num));
        std::vector<Eigen::MatrixXd> Err12(coeff1.size(), Eigen::MatrixXd::Zero(6,num));
        std::vector<Eigen::MatrixXd> Err22(coeff1.size(), Eigen::MatrixXd::Zero(6,num));
}

