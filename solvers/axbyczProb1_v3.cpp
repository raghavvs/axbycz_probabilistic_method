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

This function implements the Prob1 method in the paper
 Prerequisites on the input
   A1 is constant with B1 and C1 free
   C2 is constant with A2 adn B2 free

Inputs:
    A1, B1, C1, A2, B2, C2 : 4 x 4 x n Matrices
    opt : bool input for batchSolveXY : true - known nstd of data noise
    nstd1, nstd2 - standard deviation
Outputs:
    X_final, Y_final, Z_final - pass by reference to the solver
*/

#include <iostream>
#include <Eigen/Dense>
#include "batchSolveXY.h"
#include "rotError.h"
#include "tranError.h"

void axbyczProb1(const Eigen::Matrix4d& A1,
                 const Eigen::Matrix4d& B1,
                 const Eigen::Matrix4d& C1,
                 const Eigen::Matrix4d& A2,
                 const Eigen::Matrix4d& B2,
                 const Eigen::Matrix4d& C2,
                 int s_X,
                 int s_Y,
                 int s_Z,
                 double nstd1,
                 double nstd2,
                 bool opt,
                 Eigen::Matrix4d& X_final,
                 Eigen::Matrix4d& Y_final,
                 Eigen::Matrix4d& Z_final) {
    // Define some variables
    std::vector<Eigen::Matrix4d> X(s_X);
    //std::vector<Eigen::Matrix4d> Y(s_Y);
    std::vector<Eigen::Matrix4d> Z(s_Z);

    // Define mean and covariance matrices
    Eigen::MatrixXd MeanA1[2], MeanB1[2], MeanC1[2];
    Eigen::MatrixXd SigA1[2], SigB1[2], SigC1[2];
    Eigen::MatrixXd MeanA2[2], MeanB2[2], MeanC2[2];
    //Eigen::MatrixXd SigA2[2], SigB2[2], SigC2[2];

    // Solve for X and Z
    batchSolveXY(A1, B1, s_X, opt, nstd1, nstd2, X_final,
                 Z_final, MeanA1[0], MeanB1[0], SigA1[0],SigB1[0]);

    // Keep the candidates of Z that are SE3
    // Normally, there will be four Z \in SE3
    //std::vector<Eigen::Matrix4d> Z;
    std::vector<Eigen::Matrix4d> Z_g;
    for (const auto& z : Z_g) {
        if (z.determinant() > 0) {
            Z.push_back(z);
        }
    }

    //size_t s_Z = Z.size();

    // Solve for X
    // C2 fixed, A2 and B2 free

    // Calculate B2^-1
    //const size_t Num = sizeof(A2) / sizeof(A2[0]);
    /*size_t Num = A2.size();
    std::vector<Eigen::MatrixXd> A2_inv(Num);
    std::vector<Eigen::MatrixXd> B2_inv(Num);
    for (size_t i = 0; i < Num; ++i) {
        A2_inv[i] = A2[i].inverse();
        B2_inv[i] = B2[i].inverse();
    }*/

    size_t dim1 = 4;
    size_t dim2 = 4;
    size_t Num = A2.cols() / dim2;

    std::vector<Eigen::Matrix4d> A2_inv(Num);
    std::vector<Eigen::Matrix4d> B2_inv(Num);

    for (size_t k = 0; k < Num; ++k) {
        Eigen::Matrix4d A2_mat = A2.block(0, k*dim2, dim1, dim2);
        Eigen::Matrix4d B2_mat = B2.block(0, k*dim2, dim1, dim2);

        A2_inv[k] = A2_mat.inverse();
        B2_inv[k] = B2_mat.inverse();
    }



    // Calculate X_g : all guesses of X
    std::vector<Eigen::Matrix4d> X_g;
    batchSolveXY(A1, B1, s_X, opt, nstd1, nstd2, X_final, Z_final, MeanA1[0], MeanB1[0], SigA1[0],
                 SigB1[0]);

    // Calculate MeanB2 for computing Y later
    batchSolveXY(A1, B1, s_X, opt, nstd1, nstd2, X_final, Z_final, MeanA1[0], MeanB1[0], SigA1[0],
                 SigB1[0]);

    // Keep the candidates of X that are SE3
    // Normally, there will be four X \in SE3
    for (const auto& x : X_g) {
        if (x.determinant() > 0) {
            X.push_back(x);
        }
    }

    //size_t s_X = X.size();

    // Solve for Y
    // Compute Y using the mean equations
    /*size_t dim1 = 4;
    size_t dim2 = 4;*/
    size_t dim3 = 2 * s_X * s_Z;

    std::vector<Eigen::MatrixXd> Y(dim3, Eigen::MatrixXd::Zero(dim1, dim2));

    //Y = Zero(2 * s_X * s_Z);
    for (size_t i = 0; i < s_X; ++i) {
        for (size_t j = 0; j < s_Z; ++j) {
            std::cout << "A1: " << A1.rows() << "x" << A1.cols() << std::endl;
            std::cout << "X[i]: " << X[i].rows() << "x" << X[i].cols() << std::endl;
            std::cout << "MeanB1[0]: " << MeanB1[0].rows() << "x" << MeanB1[0].cols() << std::endl;
            std::cout << "Z[j].inverse(): " << Z[j].inverse().rows() << "x" << Z[j].inverse().cols() << std::endl;
            std::cout << "MeanC1[0].inverse(): " << MeanC1[0].inverse().rows() << "x" << MeanC1[0].inverse().cols() << std::endl;
            Y[i * s_Z + j] = (A1 * X[i] * MeanB1[0] * Z[j].inverse()) * MeanC1[0].inverse();
            Y[i * s_Z + j + s_X * s_Z] = ((MeanA2[0] * X[i] * MeanB2[0]) * Z[j].inverse()) * C2.inverse();
        }
    }

    //size_t s_Y = Y.size();
    // Find out the optimal (X, Y, Z) that minimizes cost
    Eigen::MatrixXd cost = Eigen::MatrixXd::Zero(s_X, s_Y * s_Z);
    double weight = 1.5; // weight on the translational error of the cost function
    for (size_t i = 0; i < s_X; ++i) {
        for (size_t j = 0; j < s_Z; ++j) {
            for (size_t m = 0; m < s_Y; ++m) {
                Eigen::Matrix4d left1 = A1 * X[i] * MeanB1[0];
                Eigen::Matrix4d right1 = Y[m] * MeanC1[0] * Z[j];
                double diff1 = rotError(left1, right1) + weight * tranError(left1, right1);

                Eigen::Matrix4d left2 = MeanA2[0] * X[i] * MeanB2[0];
                Eigen::Matrix4d right2 = Y[m] * C2 * Z[j];
                double diff2 = rotError(left2, right2) + weight * tranError(left2, right2);

                cost(i, j*s_Y + m) = diff1 + diff2;
            }
        }
    }

    // recover the X,Y,Z that minimizes cost
    size_t I_row;
    size_t I_col;
    double min_cost;
    min_cost = cost.minCoeff(&I_row,&I_col);

    X_final=X[I_row];

    size_t index_Z=I_col/s_Y;

    Z_final=Z[index_Z];

    size_t index_Y=I_col%s_Y;

    Y_final=Y[index_Y];
}

int main() {
    // Define some variables
    Eigen::Matrix4d A1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d B1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d C1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d A2 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d B2 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d C2 = Eigen::Matrix4d::Random();

    int s_X = 10;
    int s_Y = 10;
    int s_Z = 10;

    double nstd1 = 0.5;
    double nstd2 = 0.5;

    bool opt = true;

    // Define output variables
    Eigen::Matrix4d X_final;
    Eigen::Matrix4d Y_final;
    Eigen::Matrix4d Z_final;

    // Call the function
    axbyczProb1(A1, B1, C1, A2, B2, C2, s_X, s_Y, s_Z, nstd1, nstd2, opt,
                X_final, Y_final, Z_final);

    // Print the results
    std::cout << "X_final:\n" << X_final << std::endl;
    std::cout << "Y_final:\n" << Y_final << std::endl;
    std::cout << "Z_final:\n" << Z_final << std::endl;

    return 0;
}

/*
int main(){
    std::cout << "Build successful" << std::endl;
}*/

/*
A1: 4x4
X[i]: 4x4
MeanB1[0]: 4x4
Z[j].inverse(): 4x4
MeanC1[0].inverse(): 0x0*/