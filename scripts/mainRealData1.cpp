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

    int init_guess = 1;
    Eigen::Matrix4d X_init, Y_init, Z_init;
    Eigen::Matrix4d X_cal1, Y_cal1, Z_cal1, X_cal3, Y_cal3, Z_cal3;

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

    double err1, err3;

    std::cout << "Probability Method 1..." << std::endl;
    axbyczProb1(A1, B1, C1,
                A2, B2, C2,
                0, 0, 0,
                X_cal1, Y_cal1, Z_cal1);

    // Initial guess for iterative refinement as the results from prob 1
    if (init_guess == 3) {
        X_init = X_cal1;
        Y_init = Y_cal1;
        Z_init = Z_cal1;
    }

    // Iterative Refinement
    std::cout << "Iterative Refinement..." << std::endl;
    int num = 1;
    axbyczProb3(A1, B1, C1,
                A2, B2, C2,
                X_init, Y_init, Z_init,
                X_cal3, Y_cal3, Z_cal3,
                num);

    // Verification
    // Prob 1
    err1 = metric(A1, B1, C1, X_cal1, Y_cal1, Z_cal1) +
                        metric(A2, B2, C2, X_cal1, Y_cal1, Z_cal1);

    // Iterative refinement
    err3 = metric(A1, B1, C1, X_cal3, Y_cal3, Z_cal3) +
               metric(A2, B2, C2, X_cal3, Y_cal3,Z_cal3);


    std::ofstream outFile("results/XYZ.txt", std::ios_base::app);

    outFile << "Probability Method 1" << std::endl;
    outFile << "X_cal1: " << std::endl << X_cal1 << std::endl;
    outFile << "Y_cal1: " << std::endl << Y_cal1 << std::endl;
    outFile << "Z_cal1: " << std::endl << Z_cal1 << std::endl;
    outFile << "Probability Method 3 - Iterative Refinement" << std::endl;
    outFile << "X_cal3: " << std::endl << X_cal3 << std::endl;
    outFile << "Y_cal3: " << std::endl << Y_cal3 << std::endl;
    outFile << "Z_cal3: " << std::endl << Z_cal3 << std::endl;

    outFile.close();

    std::cout << "Error 1: " << err1 << std::endl;
    std::cout << "Error 3: " << err3 << std::endl;
}

/*
 * Error:
 *
 * Probability Method 1...
Signal: SIGSEGV (Segmentation fault)

 [mainRealData1] _mm_load_pd emmintrin.h:124
[mainRealData1] Eigen::internal::pload<double __vector(2)>(Eigen::internal::unpacket_traits<double __vector(2)>::type const*) PacketMath.h:307
[mainRealData1] Eigen::internal::ploadt<double __vector(2), 16>(Eigen::internal::unpacket_traits<double __vector(2)>::type const*) GenericPacketMath.h:463
[mainRealData1] Eigen::internal::evaluator<Eigen::PlainObjectBase<Eigen::Matrix<double, 4, 4, 0, 4, 4> > >::packet<16, double __vector(2)>(long, long) const CoreEvaluators.h:197
[mainRealData1] Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Matrix<double, 4, 4, 0, 4, 4> >, Eigen::internal::evaluator<Eigen::Matrix<double, 4, 4, 0, 4, 4> >, Eigen::internal::assign_op<double, double>, 0>::assignPacket<16, 16, double __vector(2)>(long, long) AssignEvaluator.h:652
[mainRealData1] Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Matrix<double, 4, 4, 0, 4, 4> >, Eigen::internal::evaluator<Eigen::Matrix<double, 4, 4, 0, 4, 4> >, Eigen::internal::assign_op<double, double>, 0>::assignPacketByOuterInner<16, 16, double __vector(2)>(long, long) AssignEvaluator.h:666
[mainRealData1] Eigen::internal::copy_using_evaluator_innervec_CompleteUnrolling::run AssignEvaluator.h:274
[mainRealData1] Eigen::internal::dense_assignment_loop::run AssignEvaluator.h:468
[mainRealData1] Eigen::internal::call_dense_assignment_loop<…> AssignEvaluator.h:741
[mainRealData1] Eigen::internal::Assignment::run AssignEvaluator.h:879
[mainRealData1] Eigen::internal::call_assignment_no_alias<…> AssignEvaluator.h:836
[mainRealData1] Eigen::internal::evaluator_assume_aliasing<Eigen::Matrix<double, 4, 4, 0, 4, 4>, Eigen::internal::evaluator_traits<Eigen::Matrix<double, 4, 4, 0, 4, 4> >::Shape>::value, void*>::type) AssignEvaluator.h:804
[mainRealData1] Eigen::internal::call_assignment<…> AssignEvaluator.h:782
[mainRealData1] Eigen::PlainObjectBase::_set<…> PlainObjectBase.h:714
[mainRealData1] Eigen::Matrix::operator= Matrix.h:208
[mainRealData1] axbyczProb1 axbyczProb1.h:152
[mainRealData1] main mainRealData1.cpp:101
[libc.so.6] __libc_start_main 0x00007ffff7513083
[mainRealData1] _start 0x0000555555557a5e

 */