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
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

int main() {
    // Load data
    // load('real_data/transform_ABC_unified.mat'); // load the experiment data
    // load('real_data/transform_ABC_unified_fixA.mat');
    // load('real_data/transform_ABC_unified_fixC.mat');

    // Generate Data
    // Initial guesses:
    // 1 for identity; 2(or 3) for approximate measurement from kinematics
    // data of the robot; 3 for results from Prob 1.

    int init_guess = 1;

<<<<<<< HEAD
    Matrix4d X_init;
    Matrix4d Y_init;
    Matrix4d Z_init;
    // Initial guess as Identity if (init_guess == 1) { X_init =

=======
    Eigen::Matrix4d X_init, Y_init, Z_init;

    if (init_guess == 1) {
        X_init = Eigen::Matrix4d::Identity();
        Y_init = Eigen::Matrix4d::Identity();
        Z_init = Eigen::Matrix4d::Identity();
    } else {
        // Initial guess as approx. measurement
        X_init << rotx(-M_PI / 2) * roty(M_PI / 2) * rotx(1.2 * M_PI / 180), Eigen::Vector3d(58.71,0,63.64)/1000,
                0,0,0,1;
        Y_init << rotz(M_PI)*rotz(M_PI/4), Eigen::Vector3d(400,0,0)/1000,
                0,0,0,1;
        Z_init << rotz(M_PI), Eigen::Vector3d(0,0,10)/1000,
                0,0,0,1;
    }

    // Choice of scramble rate
    bool isRandPerm = true;
    std::vector<int> r = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    for (int rk = 0; rk < r.size(); rk++) {
        // convert cells to matrices
        std::vector<Eigen::MatrixXd> A;
        std::vector<Eigen::MatrixXd> B;
        std::vector<Eigen::MatrixXd> C;
        std::vector<Eigen::MatrixXd> Bp;

        for (int i = 0; i < headTf1.size(); i++) {
            // Inputs for Iterative Refinement
            Eigen::MatrixXd A1_i;
            Eigen::MatrixXd B1_i;
            Eigen::MatrixXd C1_i;
            std::tie(A1_i,B1_i,C1_i) = convertCell2Mat(headTf1[i], handTf1[i], tagTf1[i]);

            Eigen::MatrixXd A2_i;
            Eigen::MatrixXd B2_i;
            Eigen::MatrixXd C2_i;
            std::tie(A2_i,B2_i,C2_i) = convertCell2Mat(headTf2[i], handTf2[i], tagTf2[i]);

            if (isRandPerm) {
                Bp1[i] = scrambleData(B1[i], r[rk]);
                Bp2[i] = scrambleData(B2[i], r[rk]);
            }

            // Inputs for Wang
            Eigen::MatrixXd AA = Eigen::MatrixXd(A1_i.rows(), A1_i.cols(), 2);
            AA << A1_i, A2_i;
            Eigen::MatrixXd BB = Eigen::MatrixXd(B1_i.rows(), B1_i.cols(), 2);
            BB << B1_i, B2_i;
            Eigen::MatrixXd CC = Eigen::MatrixXd(C1_i.rows(), C1_i.cols(), 2);
            CC << C1_i, C2_i;

            A.push_back(AA);
            B.push_back(BB);
            C.push_back(CC);
        }

        // Inputs for Prob 1
        Eigen::MatrixXd AA1;
        Eigen::MatrixXd BB1;
        Eigen::MatrixXd CC1;
        std::tie(AA1,BB1,CC1) = convertCell2Mat(headTf1[0], handTf1[0], tagTf1[0]);

        Eigen::MatrixXd AA2;
        Eigen::MatrixXd BB2;
        Eigen::MatrixXd CC2;

        std::tie(AA2,BB2,CC2) = convertCell2Mat(headTf2[0], handTf2[0], tagTf2[0]);
        if (isRandPerm) {
            Eigen::MatrixXd BBp1 = scrambleData(BB1, r[rk]);
            Eigen::MatrixXd BBp2 = scrambleData(BB2, r[rk]);
            Eigen::MatrixXd Bp = scrambleData(B, r[rk]);
        }

        // Prob 1
        std::cout << "Probabilistic Method 1..." << std::endl;
        Eigen::Matrix4d X_cal1;
        Eigen::Matrix4d Y_cal1;
        Eigen::Matrix4d Z_cal1;
        std::tie(X_cal1,Y_cal1,Z_cal1) = axbyczProb1(AA1.slice(0), BBp1, CC1,
                                                     AA2, BBp2, CC2.slice(0), 0, 0, 0);

        // Initial guess for iterative refinement as the results from prob 1
        if (init_guess == 3) {
            X_init = X_cal1;
            Y_init = Y_call;
            Z_init = Z_call;
        }

        // Iterative Refinement
        std::cout << "Iterative Refinement..." << std::endl;
        Eigen::Matrix4d X_cal2;
        Eigen::Matrix4d Y_cal2;
        Eigen::Matrix4d Z_cal2;
        int num2_rk;
        std::tie(X_cal2,Y_cal2,Z_cal2,num2_rk) = axbyczProb3(A1, Bp1, C1,
                                                             A2, Bp2, C2, X_init, Y_init, Z_init);
        num2[rk] = num2_rk;

        // Call traditional AXB=YCZ algorithm to solve for X, Y and Z given A, B and C
        std::cout << "Wang Method..." << std::endl;
        Eigen::Matrix4d X_cal3;
        Eigen::Matrix4d Y_cal3;
        Eigen::Matrix4d Z_cal3;
        int num3_rk;
        std::tie(X_cal3,Y_cal3,Z_cal3,num3_rk) = Wang_AXBYCZ(A, Bp, C,
                                                             X_init, Y_init, Z_init);
        num3[rk] = num3_rk;

        // Verification
        // Prob 1
        err1[rk] = metric(A1,B1,C1,X_cal1,Y_cal1,Z_cal1) +
                   metric(A2,B2,C2,X_cal1,Y_cal1,Z_cal1);

        // Iterative refinement
        err2[rk] = metric(A1,B1,C1,X_cal2,Y_cal2,Z_cal2) +
                   metric(A2,B2,C2,X_cal2,Y_cal2,Z_cal2);

        // Wang method
        err3[rk] = metric(A1,B1,C1,X_cal3,Y_cal3,Z_cal3) +
                   metric(A2,B2,C2,X_cal3,Y_cal3,Z_cal3);
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

    //-- Supporting functions --//

    // Convert data to 3d matrices
    std::tuple<Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXd>
            convertCell2Mat(std::vector<Eigen::Matrix4d> headTf,
                            std::vector<Eigen::Matrix4d> handTf,
                            std::vector<Eigen::Matrix4d> tagTf) {
        Eigen::MatrixXd A(headTf[0].rows(), headTf[0].cols(), headTf.size());
        Eigen::MatrixXd B(tagTf[0].rows(), tagTf[0].cols(), tagTf.size());
        Eigen::MatrixXd C(handTf[0].rows(), handTf[0].cols(), handTf.size());
        for (int i = 0; i < headTf.size(); ++i) {
            A.slice(i) = headTf[i];
            B.slice(i) = tagTf[i];
            C.slice(i) = handTf[i];
        }
        return std::make_tuple(A,B,C);
    }

    // Metric for computing errors
    double metric(std::vector<Eigen::Matrix4d> A,
                  std::vector<Eigen::Matrix4d> B,
                  std::vector<Eigen::Matrix4d> C,
                  Eigen::Matrix4d X,
                  Eigen::Matrix4d Y,
                  Eigen::Matrix4d Z) {
        double diff = 0;
        int N = 0;
        for (int i = 0; i < A.size(); ++i) {
            for (int j = 0; j < A[i].size(); ++j) {
                diff += (A[i].slice(j)*X*B[i].slice(j)-Y*C[i].slice(j)*Z).norm();
                N++;
            }
        }
        return diff/N;
    }
>>>>>>> dev_mac
}