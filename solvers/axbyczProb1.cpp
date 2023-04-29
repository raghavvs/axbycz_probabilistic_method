/*
DESCRIPTION:

The provided code implements a set of functions that solve a robotics problem
involving the transformation matrices of multiple coordinate frames. Specifically,
the functions solve for the transformations between three coordinate frames
(A, B, and C) given the transformations between A and B and between A and C.
The functions use the Eigen library to perform matrix operations such as inversion
and SVD decomposition. The main function (axbyczProb1) calls the other two functions
(batchSolveXY and randSE3) to generate a set of random transformations and iteratively
select those that satisfy certain constraints, in order to estimate the desired transformations.

Input:
A1, B1, C1, A2, B2, C2: Matrices - dim 4x4
opt: bool
nstd1, nst2: standard deviation
Output:
X_final, Y_final, Z_final: Matrices - dim 4x4

In the case of two robotic arms:
 A - robot 1's base to end effector transformation (forward kinematics)
 B - camera to calibration target transformation
 C - robot 2's base to end effector transformation (forward kinematics)
 X - end effector of robot 1 to camera transformation
 Y - robot 1's base to robot 2's base transformation
 Z - end effector of robot 2 to calibration target transformation
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <eigen3/Eigen/Dense>
#include "batchSolveXY.h"
#include "loadMatrices.h"
#include "rotError.h"
#include "tranError.h"

void axbyczProb1(const std::vector<Eigen::Matrix4d>& A1,
                 const std::vector<Eigen::Matrix4d>& B1,
                 const std::vector<Eigen::Matrix4d>& C1,
                 const std::vector<Eigen::Matrix4d>& A2,
                 const std::vector<Eigen::Matrix4d>& B2,
                 const std::vector<Eigen::Matrix4d>& C2,
                 bool opt,
                 double nstd1,
                 double nstd2,
                 Eigen::Matrix4d& X_final,
                 Eigen::Matrix4d& Y_final,
                 Eigen::Matrix4d& Z_final){

    //// A1 is constant with B1 and C1 free
    Eigen::Matrix4d A1_fixed = A1[0];

    //// C2 is constant with A2 and B2 free
    Eigen::Matrix4d C2_fixed = C2[0];

    //// Solve for Z
    std::vector<Eigen::Matrix4d> Z_g, X_dummy, Y_dummy;
    Eigen::Matrix4d MeanC1, MeanB1, MeanA2, MeanB2, MeanA2_inv, MeanB2_inv;
    Eigen::Matrix<double, 6, 6> SigC1, SigB1, SigA2, SigB2;

    // Calculate MeanC1
    batchSolveXY(C1, B1, opt, nstd1, nstd2, Z_g, Y_dummy,
                 MeanC1, MeanB1, SigC1, SigB1);

    std::cout << "Z_g[last]: " << std::endl << Z_g[1] << std::endl;

    std::vector<Eigen::Matrix4d> Z;
    for (const auto& z : Z_g) {
        if (z.determinant() > 0) {
            Z.push_back(z);
        }
    }

    int s_Z = Z.size();

    // all good

    // Calculate B2_inv
    int Num = A2.size();
    std::vector<Eigen::Matrix4d> A2_inv(Num), B2_inv(Num);
    for (int i = 0; i < Num; ++i) {
        A2_inv[i] = A2[i].inverse();
        B2_inv[i] = B2[i].inverse();
    }

    //// Solve for X
    std::vector<Eigen::Matrix4d> X_g;

    // Calculate MeanA2
    batchSolveXY(A2, B2_inv, opt, nstd1, nstd2, X_g, Y_dummy,
                 MeanA2, MeanB2_inv, SigA2, SigB2);

    std::vector<Eigen::Matrix4d> X;
    for (const auto& x : X_g) {
        if (x.determinant() > 0) {
            X.push_back(x);
        }
    }

    size_t s_X = X.size();

    // Calculate MeanB2 for computing Y later
    batchSolveXY(A2_inv, B2, opt, nstd1, nstd2, X_dummy, Y_dummy,
                 MeanA2_inv, MeanB2, SigA2, SigB2);

    // all good

    // Solve for Y
    std::vector<Eigen::Matrix4d> Y(2 * s_X * s_Z);
    for (int i = 0; i < s_X; ++i) {
        for (int j = 0; j < s_Z; ++j) {
            // ignore element wise multiplication form in MATLAB version
            Y[(i * s_Z) + j] = (A1_fixed * X[i] * MeanB1 * Z[j].inverse()) * MeanC1.inverse();
            Y[(i * s_Z) + j + s_X * s_Z] = (MeanA2 * X[i] * MeanB2 * Z[j].inverse()) * C2_fixed.inverse();
        }
    }

    int s_Y = Y.size();

    std::cout << s_X << std::endl;
    std::cout << s_Y << std::endl;
    std::cout << s_Z << std::endl;

    std::cout << "X[last]: " << std::endl << X[s_X-1] << std::endl;
    std::cout << "Y[last]: " << std::endl << Y[s_Y-1] << std::endl;
    std::cout << "Z[last]: " << std::endl << Z[s_Z-1] << std::endl;

    // all good

    // Find the optimal (X, Y, Z) that minimizes cost
    Eigen::MatrixXd cost(s_X, s_Y * s_Z);
    double weight = 1.5;
    double min_cost = std::numeric_limits<double>::max();
    int min_i = 0, min_j = 0, min_m = 0;

    for (size_t i = 0; i < s_X; ++i) {
        for (size_t j = 0; j < s_Z; ++j) {
            for (size_t m = 0; m < s_Y; ++m) {
                Eigen::Matrix4d left1 = A1_fixed * X[i] * MeanB1;
                Eigen::Matrix4d right1 = Y[m] * MeanC1 * Z[j];

                double diff1 = rotError(left1, right1) + weight *
                                                         tranError(left1, right1);

                Eigen::Matrix4d left2 = MeanA2 * X[i] * MeanB2;
                Eigen::Matrix4d right2 = Y[m] * C2_fixed * Z[j];

                double diff2 = rotError(left2, right2) + weight *
                                                         tranError(left2, right2);

                double current_cost = std::abs(diff1) + std::abs(diff2);

                cost(i, j * s_Y + m) = current_cost;
                if (current_cost < min_cost) {
                    min_cost = current_cost;
                    min_i = i;
                    min_j = j;
                    min_m = m;
                }
            }
        }
    }

    //// Recover the X, Y, Z that minimize cost
    X_final = X[min_i];
    Y_final = Y[min_m];
    Z_final = Z[min_j];

    std::cout << min_i << std::endl;
    std::cout << min_j << std::endl;
    std::cout << min_m << std::endl;
}

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

    // Set opt, nstd1, and nstd2
    bool opt = true;
    double nstd1 = 0.01;
    double nstd2 = 0.01;

    // Resize the vectors to keep only the first 10 matrices
    int new_size = 10;
    A1.resize(new_size);
    B1.resize(new_size);
    C1.resize(new_size);
    A2.resize(new_size);
    B2.resize(new_size);
    C2.resize(new_size);

    // Call the axbyczProb1 function
    Eigen::Matrix4d X_final, Y_final, Z_final;
    axbyczProb1(A1, B1, C1, A2, B2, C2, opt, nstd1, nstd2, X_final, Y_final, Z_final);

    // Display results
    std::cout << "X_final:\n" << X_final << std::endl;
    std::cout << "Y_final:\n" << Y_final << std::endl;
    std::cout << "Z_final:\n" << Z_final << std::endl;

    return 0;
}

/*
int main() {
    // Create deterministic input matrices
    int num_matrices = 2;
    std::vector<Eigen::Matrix4d> A1(num_matrices), B1(num_matrices), C1(num_matrices);
    std::vector<Eigen::Matrix4d> A2(num_matrices), B2(num_matrices), C2(num_matrices);

    // Fill in the input matrices with specific examples
    A1[0] << 1, 2, 3, 1,
            0, 1, 0, 2,
            0, 0, 1, 3,
            0, 0, 0, 1;
    A1[1] << 1, 1, 3, 2,
            0, 1, 0, 3,
            0, 0, 1, 1,
            0, 0, 0, 1;

    B1 = A1;
    C1 = A1;
    A2 = A1;
    B2 = A1;
    C2 = A1;

    // Set opt, nstd1, and nstd2
    bool opt = true;
    double nstd1 = 0.01;
    double nstd2 = 0.01;

    // Call the axbyczProb1 function
    Eigen::Matrix4d X_final, Y_final, Z_final;
    axbyczProb1(A1, B1, C1, A2, B2, C2, opt, nstd1, nstd2, X_final, Y_final, Z_final);

    // Display results
    std::cout << "X_final:\n" << X_final << std::endl;
    std::cout << "Y_final:\n" << Y_final << std::endl;
    std::cout << "Z_final:\n" << Z_final << std::endl;

    return 0;
}*/
