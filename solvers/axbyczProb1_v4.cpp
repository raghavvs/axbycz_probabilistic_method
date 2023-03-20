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

void axbyczProb1(const Eigen::Matrix4d& A1,
                 const Eigen::Matrix4d& B1,
                 const Eigen::Matrix4d& C1,
                 const Eigen::Matrix4d& A2,
                 const Eigen::Matrix4d& B2,
                 const Eigen::Matrix4d& C2,
                 int len,
                 bool opt,
                 double nstd1,
                 double nstd2,
                 std::vector<Eigen::MatrixXd> &X_final,
                 std::vector<Eigen::MatrixXd> &Y_final,
                 std::vector<Eigen::MatrixXd> &Z_final) {
    // This function implements the Prob1 method in the paper
    // Prerequisites on the input
    //   A1 is constant with B1 and C1 free
    //   C2 is constant with A2 adn B2 free

    auto A = A1.block(0, 0, 4, 4);
    auto C = C2.block(0, 0, 4, 4);

    // ------ Solve for Z -------- //
    // A fixed, B and C free

    // ------ using probability methods ------
    // calculate Z_g : all guesses of Z
    std::vector<Eigen::MatrixXd> Z_g(len);
    Eigen::MatrixXd MeanC;
    Eigen::MatrixXd MeanB;
    Eigen::MatrixXd SigA;
    Eigen::MatrixXd SigB;

    batchSolveXY(C1,B1,len,opt,nstd1,nstd2,Z_g,Y_final,MeanC,MeanB,SigA,SigB);
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
    std::cout << "X_final: " << Y_final[0] << std::endl;

    return 0;
}
