/*
DESCRIPTION:

This function implements the Prob2 method in the paper
Prerequisites on the input
  A1 is constant with B1 and C1 free
  C2 is constant with A1 adn B1 free
  B3 is constant with A3 and C3 free

Input:
    A1, B1, C1, A2, B2, C2: Matrices - dim 4x4
    opt: bool
    nstd1, nst2: standard deviation
Output:
    X_final, Y_final, Z_final: Matrices - dim 4x4
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include "batchSolveXY.h"
#include "rotError.h"
#include "tranError.h"

void axbyczProb2(const Eigen::Matrix4d &A1,
                 const Eigen::Matrix4d &B1,
                 const Eigen::Matrix4d &C1,
                 const Eigen::Matrix4d &A2,
                 const Eigen::Matrix4d &B2,
                 const Eigen::Matrix4d &C2,
                 const Eigen::Matrix4d &A3,
                 const Eigen::Matrix4d &B3,
                 const Eigen::Matrix4d &C3,
                 std::vector<Eigen::Matrix4d> &X_final,
                 std::vector<Eigen::Matrix4d> &Y_final,
                 std::vector<Eigen::Matrix4d> &Z_final) {

    // A1 fixed, B1 and C1 free
    std::cout << "Running Prob2 method ... \n";

    int len = 8;
    double nstd1 = 0;
    double nstd2 = 0;
    bool opt = false;

    std::vector<Eigen::Matrix4d> Z_g(len);
    Eigen::MatrixXd MeanA, MeanB, MeanC, SigA, SigB, SigC;

    std::vector<Eigen::Matrix4d> A1_vec(len), B1_vec(len), C1_vec(len),
                                A2_vec(len), B2_vec(len), C2_vec(len),
                                A3_vec(len), B3_vec(len), C3_vec(len);
    for (int i = 0; i < len; ++i) {
        A1_vec[i] = A1.block(0, 0, 4, 4);
        B1_vec[i] = B1.block(0, 0, 4, 4);
        C1_vec[i] = C1.block(0, 0, 4, 4);
        A2_vec[i] = A2.block(0, 0, 4, 4);
        B2_vec[i] = B2.block(0, 0, 4, 4);
        C2_vec[i] = C2.block(0, 0, 4, 4);
        A3_vec[i] = A3.block(0, 0, 4, 4);
        B3_vec[i] = B3.block(0, 0, 4, 4);
        C3_vec[i] = C3.block(0, 0, 4, 4);
    }

    //// ------ using probability methods ------
    // calculate Z_g : all guesses of Z
    //// ------ Solve for Z -------- //
    // A1 fixed, B1 and C1 free

    Eigen::MatrixXd MeanB1, MeanC1, MeanA2;

    batchSolveXY(C1_vec, B1_vec, len, opt,nstd1,nstd2,Z_g,Y_final,
                 MeanC1,MeanB1,SigC,SigB);

    // Keep the candidates of Z that are SE3
    // Normally there will be four Z \in SE3
    int Z_index = 0;
    for (int i = 0; i < Z_g.size(); ++i) {
        if (Z_g[i].determinant() > 0) {
            Z_final.push_back(Z_g[i]);
            ++Z_index;
        }
    }

    std::cout << "works till here - solve for Z? - YES" << std::endl;

    //// ------ Solve for X -------- //
    // C2 fixed, A2 and B2 free

    // ------ Calculate B2^-1 -------
    int Num = A2.size();
    std::vector<Eigen::Matrix4d> A2_inv(Num), B2_inv(Num);
    for (int i = 0; i < Num; ++i) {
        A2_inv[i] = A2_vec[i].inverse();
        B2_inv[i] = B2_vec[i].inverse();
    }

    // ------ using probability methods ------
    // calculate X_g : all guesses of X
    std::vector<Eigen::Matrix4d> X_g(len);
    batchSolveXY(A2_vec, B2_inv, len, opt, nstd1, nstd2, X_g, Y_final,
                 MeanA2, MeanB, SigC, SigB);

    // Calculate MeanB for computing Y later
    // Note: can be further simplified by using only the distribution function
    Eigen::MatrixXd MeanB2;
    batchSolveXY(A2_inv, B2_vec, len, opt, nstd1, nstd2, X_g, Y_final,
                 MeanA, MeanB2, SigC, SigB);

    // Keep the candidates of X that are SE3å
    // Normally there will be four X \in SE3
    int X_index = 0;
    std::vector<Eigen::Matrix4d> X;
    for (int i = 0; i < X_g.size(); ++i) {
        if (X_g[i].determinant() > 0) {
            X.push_back(X_g[i]);
            ++X_index;
        }
    }

    std::cout << "works till here - solve for Z and X? - YES" << std::endl;

    //// ------ Solve for Y -------- //
    // B3 fixed, A3 and C3 free

    // ------ Calculate B2^-1 -------
    std::vector<Eigen::Matrix4d> A3_inv(len);
    std::vector<Eigen::Matrix4d> C3_inv(len);
    for (int i = 0; i < len; i++) {
        A3_inv[i] = A3_vec[i].inverse();
        C3_inv[i] = C3_vec[i].inverse();
    }

    // ------ using probability methods ------
    // calculate X_g : all guesses of X
    std::vector<Eigen::Matrix4d> Y_g_inv;
    batchSolveXY(C3_inv, A3_inv, len, opt, nstd1, nstd2, Y_g_inv, Y_final,
                 MeanC, MeanA, SigC, SigA);

    // Calculate MeanA2 and MeanC2 for the cost function later
    Eigen::MatrixXd MeanC3, MeanA3;
    batchSolveXY(C3_vec, A3_vec, len, opt, nstd1, nstd2, Y_g_inv, Y_final,
                 MeanC3, MeanA3, SigC, SigA);

    // Keep the candidates of Y that are SE3
    int Y_index = 1;
    std::vector<Eigen::Matrix4d> Y(4, Eigen::Matrix4d::Zero());
    for (int i = 0; i < Y_g_inv.size(); i++) {
        if (Y_g_inv[i].determinant() > 0) {
            Y[Y_index] = Y_g_inv[i].inverse();
            Y_index++;
        }
    }

    std::cout << "works till here - solve for Z, X and Y? - YES " << std::endl;

    //// Find out the optimal (X, Y, Z) that minimizes cost

    int s_Z = Z_final.size();
    int s_X = X.size();
    int s_Y = Y.size();


    Eigen::MatrixXd cost = Eigen::MatrixXd::Zero(s_X, s_Y * s_Z);
    double weight = 1.8; // weight on the translational error of the cost function

    for (int i = 0; i < s_X; i++) {
        for (int j = 1; j < s_Y; j++) {
            for (int p = 0; p < s_Z; p++) {
                Eigen::MatrixXd left1 = A1 * X[i] * MeanB1;
                Eigen::MatrixXd right1 = Y[j] * MeanC1 * Z_final[p];
                double diff1 = rotError(left1, right1) + weight * tranError(left1, right1);

                Eigen::MatrixXd left2 = MeanA2 * X[i] * MeanB2;
                Eigen::MatrixXd right2 = Y[j] * C2 * Z_final[p];
                double diff2 = rotError(left2, right2) + weight * tranError(left2, right2);

                Eigen::MatrixXd left3 = MeanA3 * X[i] * B3;
                Eigen::MatrixXd right3 = Y[j] * MeanC3 * Z_final[p];
                double diff3 = rotError(left3, right3) + weight * tranError(left3, right3);

                // How to better design the cost function is still an open
                // question
                cost(i, (j - 1) * s_Z + p) = std::norm(diff1) + std::norm(diff2) + std::norm(diff3);
            }
        }
    }

    std::cout << "works till here - optimal cost? - YES" << std::endl;

    std::cout << "cost.rows(): " << cost.rows() << "cost.cols(): " << cost.cols() << std::endl;

    //// recover the X,Y,Z that minimizes cost

    Eigen::Map<Eigen::VectorXd> cost_vec(cost.data(), cost.size());
    Eigen::MatrixXd::Index minIndex;
    cost_vec.minCoeff(&minIndex);
    int I1 = static_cast<int>(minIndex);
    int I_row = I1 / s_Y;
    int I_col = I1 % s_Y;

    Eigen::Matrix4d X_final_ = X[I_row];

    int index_Y;
    if (I_col % s_Y > 0) {
        if (I_col % s_Y == s_Z) {
            index_Y = s_Z;
        }
        index_Y = std::floor(I_col / s_Y) + 1;
    } else {
        index_Y = std::floor(I_col / s_Y);
    }

    Eigen::Matrix4d Y_final_ = Y[index_Y]; // final Y

    int index_Z;
    if (I_col % 4 > 0) {
        index_Z = I_col % 4;
    } else {
        index_Z = 4;
    }

    Eigen::Matrix4d Z_final_ = Z_final[index_Z]; // final Z

    std::cout << "works till here - recover X,Y,Z final? - YES" << std::endl;

    std::cout << "X_final_: " << std::endl << X_final_ << std::endl;
    std::cout << "Y_final_: " << std::endl << Y_final_ << std::endl;
    std::cout << "Z_final_: " << std::endl << Z_final_ << std::endl;
}

int main() {
    Eigen::Matrix4d A1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d B1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d C1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d A2 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d B2 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d C2 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d A3 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d B3 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d C3 = Eigen::Matrix4d::Random();

    std::vector<Eigen::Matrix4d> X_final;
    std::vector<Eigen::Matrix4d> Y_final;
    std::vector<Eigen::Matrix4d> Z_final;

    axbyczProb2(A1,B1,C1,A2,B2,C2,A3,B3,C3,X_final,Y_final,Z_final);

    std::cout << "Build successful? - YES" <<std::endl;

    return 0;
}

/*
Output:
len: 8
works till here - solve for Z? - YES
works till here - solve for Z and X? - YES
Y:
-2.58501 0.864086  1.20978 -1.80958
 2.31146 0.523348 0.238907 0.646305
-1.47282 -0.21759  2.77298 -1.66964
-5.26496 -2.34186 0.531769 -1.83458
works till here - solve for Z, X and Y? - YES
works till here - optimal cost? - YES
X_final_:
          1           0           0           0
          0   -0.894427   -0.447214 4.80553e-49
          0    0.447214   -0.894427 9.61107e-49
          0           0           0           1
Y_final_:
 -0.808833  -0.843726   -2.20955    1.53485
  0.458386 -0.0436236   0.150767  -0.366332
  -1.15638   -2.61803   -1.09865    1.13648
 -0.179155  -0.866537   0.647998    1.02978
Z_final_:
         -1           0           0 3.38334e-32
          0    0.786109   -0.618088 4.90051e-32
          0   -0.618088   -0.786109  1.0926e-31
          0           0           0           1
works till here - recover X,Y,Z final? - YES
Build successful? - YES
 */