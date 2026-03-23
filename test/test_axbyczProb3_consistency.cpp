/*
 * Consistency test for axbyczProb3 solver.
 * Generates valid SE(3) test data satisfying AXB = YCZ,
 * runs the solver multiple times, and verifies identical results.
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "axbyczProb3.h"

// Create a proper SE(3) matrix from angle-axis + translation
Eigen::Matrix4d makeSE3(double angle, const Eigen::Vector3d &axis,
                        const Eigen::Vector3d &trans) {
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.block<3,3>(0,0) = Eigen::AngleAxisd(angle, axis.normalized()).toRotationMatrix();
    T.block<3,1>(0,3) = trans;
    return T;
}

// SE(3) inverse
Eigen::Matrix4d SE3inverse(const Eigen::Matrix4d &T) {
    Eigen::Matrix4d Tinv = Eigen::Matrix4d::Identity();
    Eigen::Matrix3d R = T.block<3,3>(0,0);
    Eigen::Vector3d t = T.block<3,1>(0,3);
    Tinv.block<3,3>(0,0) = R.transpose();
    Tinv.block<3,1>(0,3) = -R.transpose() * t;
    return Tinv;
}

int main() {
    // Ground truth X, Y, Z as proper SE(3) matrices
    Eigen::Matrix4d X_true = makeSE3(0.3, Eigen::Vector3d(1, 0.5, 0.2), Eigen::Vector3d(0.1, -0.2, 0.3));
    Eigen::Matrix4d Y_true = makeSE3(0.5, Eigen::Vector3d(0.2, 1, 0.3), Eigen::Vector3d(-0.1, 0.4, 0.1));
    Eigen::Matrix4d Z_true = makeSE3(0.7, Eigen::Vector3d(0.3, 0.2, 1), Eigen::Vector3d(0.2, 0.1, -0.3));

    const int N = 20;
    std::vector<Eigen::Matrix4d> A1(N), B1(N), C1(N);
    std::vector<Eigen::Matrix4d> A2(N), B2(N), C2(N);

    // Generate consistent data: AXB = YCZ => B = (AX)^(-1) * Y * C * Z
    for (int i = 0; i < N; ++i) {
        double a1 = 0.1 * (i + 1);
        Eigen::Vector3d ax1((i+1)*0.1, (i+2)*0.05, (i+3)*0.03);
        Eigen::Vector3d t1((i+1)*0.01, -(i+2)*0.02, (i+3)*0.01);
        A1[i] = makeSE3(a1, ax1, t1);

        double c1 = 0.15 * (i + 1);
        Eigen::Vector3d cx1((i+3)*0.07, (i+1)*0.09, (i+2)*0.05);
        Eigen::Vector3d tc1(-(i+1)*0.02, (i+2)*0.01, (i+3)*0.03);
        C1[i] = makeSE3(c1, cx1, tc1);

        // B = (AX)^{-1} * Y * C * Z
        B1[i] = SE3inverse(A1[i] * X_true) * Y_true * C1[i] * Z_true;

        // Second set with different data
        double a2 = 0.12 * (i + 1);
        Eigen::Vector3d ax2((i+2)*0.08, (i+3)*0.06, (i+1)*0.04);
        Eigen::Vector3d t2((i+2)*0.015, (i+1)*0.025, -(i+3)*0.01);
        A2[i] = makeSE3(a2, ax2, t2);

        double c2 = 0.18 * (i + 1);
        Eigen::Vector3d cx2((i+1)*0.06, (i+3)*0.04, (i+2)*0.08);
        Eigen::Vector3d tc2((i+3)*0.02, -(i+1)*0.03, (i+2)*0.015);
        C2[i] = makeSE3(c2, cx2, tc2);

        B2[i] = SE3inverse(A2[i] * X_true) * Y_true * C2[i] * Z_true;
    }

    // Use initial guesses close to ground truth (simulating axbyczProb1 output)
    // axbyczProb3 is an iterative refinement solver, needs good starting point
    Eigen::Matrix4d Xinit_base = X_true;
    Eigen::Matrix4d Yinit_base = Y_true;
    Eigen::Matrix4d Zinit_base = Z_true;

    // Run solver 3 times and store results
    std::vector<Eigen::Matrix4d> X_results(3), Y_results(3), Z_results(3);

    for (int run = 0; run < 3; ++run) {
        Eigen::Matrix4d Xinit = Xinit_base;
        Eigen::Matrix4d Yinit = Yinit_base;
        Eigen::Matrix4d Zinit = Zinit_base;
        int num = 0;

        axbyczProb3(A1, B1, C1, A2, B2, C2,
                    Xinit, Yinit, Zinit,
                    X_results[run], Y_results[run], Z_results[run], num);

        std::cout << "=== Run " << (run + 1) << " (iterations: " << num << ") ===" << std::endl;
        std::cout << "X_cal:\n" << X_results[run] << "\n" << std::endl;
        std::cout << "Y_cal:\n" << Y_results[run] << "\n" << std::endl;
        std::cout << "Z_cal:\n" << Z_results[run] << "\n" << std::endl;
    }

    // Check consistency
    bool consistent = true;
    for (int run = 1; run < 3; ++run) {
        double x_diff = (X_results[run] - X_results[0]).norm();
        double y_diff = (Y_results[run] - Y_results[0]).norm();
        double z_diff = (Z_results[run] - Z_results[0]).norm();

        if (x_diff > 1e-10 || y_diff > 1e-10 || z_diff > 1e-10) {
            std::cout << "INCONSISTENCY detected in run " << (run + 1) << ":" << std::endl;
            std::cout << "  X diff: " << x_diff << std::endl;
            std::cout << "  Y diff: " << y_diff << std::endl;
            std::cout << "  Z diff: " << z_diff << std::endl;
            consistent = false;
        }
    }

    // Check accuracy vs ground truth
    double x_err = (X_results[0] - X_true).norm();
    double y_err = (Y_results[0] - Y_true).norm();
    double z_err = (Z_results[0] - Z_true).norm();
    std::cout << "\nAccuracy vs ground truth:" << std::endl;
    std::cout << "  X error: " << x_err << std::endl;
    std::cout << "  Y error: " << y_err << std::endl;
    std::cout << "  Z error: " << z_err << std::endl;

    if (consistent) {
        std::cout << "\nRESULT: All 3 runs produced IDENTICAL results - PASS" << std::endl;
    } else {
        std::cout << "\nRESULT: Runs produced DIFFERENT results - FAIL" << std::endl;
    }

    return consistent ? 0 : 1;
}
