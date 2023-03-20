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
*/

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "batchSolveXY.h"

void axbyczProb1(const Eigen::Matrix4d &A1,
                 const Eigen::Matrix4d &B1,
                 const Eigen::Matrix4d &C1,
                 const Eigen::Matrix4d &A2,
                 const Eigen::Matrix4d &B2,
                 const Eigen::Matrix4d &C2,
                 int len,
                 bool opt,
                 double nstd1,
                 double nstd2,
                 std::vector<Eigen::MatrixXd> &X_final,
                 std::vector<Eigen::MatrixXd> &Y_final,
                 std::vector<Eigen::MatrixXd> &Z_final) {

    //   A1 is constant with B1 and C1 free
    //   C2 is constant with A2 and B2 free

    auto A = A1.block(0, 0, 4, 4);
    auto C = C2.block(0, 0, 4, 4);

    //// ------ Solve for Z -------- //
    // A1 fixed, B1 and C1 free

    //// ------ using probability methods ------
    // calculate Z_g : all guesses of Z
    std::vector<Eigen::MatrixXd> Z_g(len);
    Eigen::MatrixXd MeanC, MeanB, SigA, SigB;

    batchSolveXY(C1,B1, len, opt,nstd1,nstd2,Z_g,Y_final,
                 MeanC,MeanB,SigA,SigB);

    // Keep the candidates of Z that are SE3
    // Normally there will be four Z \in SE3
    int Z_index = 0;
    for (int i = 0; i < len; i++) {
        if (Z_g[i].determinant() > 0) {
            Z_final.push_back(Z_g[i]);
            Z_index++;
        }
    }

    int s_Z = Z_final.size();

    std::cout << "works till here?" << std::endl;

    // Solve for X
    // C2 fixed, A2 and B2 free

    // Calculate B2^-1
    size_t dim1 = 4;
    size_t dim2 = 4;
    int Num = A2.size();
    std::vector<Eigen::MatrixXd> A2_inv(Num), B2_inv(Num);
    for (size_t k = 0; k < Num; ++k) {
        Eigen::Matrix4d A2_mat = A2.block(0, k*dim2, dim1, dim2);
        Eigen::Matrix4d B2_mat = B2.block(0, k*dim2, dim1, dim2);

        A2_inv[k] = A2_mat.inverse();
        B2_inv[k] = B2_mat.inverse();
    }

    // using probability methods
    // calculate X_g : all guesses of X
    std::vector<Eigen::MatrixXd> X_g(len);
    Eigen::MatrixXd MeanA2;
    batchSolveXY(A2, B2_inv, len, opt, nstd1, nstd2, X_g, Y_final,
                 MeanA2, MeanB, SigA, SigB);

    // Calculate MeanB for computing Y later
    batchSolveXY(A2_inv,B2,opt,nstd1,nstd2,X_g,Y_final,
                 MeanA2,MeanB,SigA,SigB);

    // Keep the candidates of X that are SE3
    // Normally there will be four X \in SE3
    int X_index = 0;
    for (int i = 0; i < Num; i++) {
        if (X_g[i].determinant() > 0) {
            X_final.push_back(X_g[i]);
            X_index++;
        }
    }

    int s_X = X_final.size();

}

int main() {
    Eigen::Matrix4d A1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d B1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d C1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d A2 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d B2 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d C2 = Eigen::Matrix4d::Random();

    int len = 2;
    bool opt = true;
    double nstd1 = 0.5;
    double nstd2 = 0.5;

    std::vector<Eigen::MatrixXd> X_final(len);
    std::vector<Eigen::MatrixXd> Y_final(len);
    std::vector<Eigen::MatrixXd> Z_final(len);

    axbyczProb1(A1,B1,C1,A2,B2,C2,len,opt,nstd1,nstd2,X_final,Y_final,Z_final);

    std::cout << "Build successful" <<std::endl;
    std::cout << "Y_final: " << std::endl << Y_final[0] << std::endl;

    return 0;
}

