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

using namespace Eigen;

void axbyczProb1(MatrixXd A1, MatrixXd B1, MatrixXd C1, MatrixXd A2, MatrixXd B2, MatrixXd C2, std::string opt, double nstd1, double nstd2, MatrixXd& X_final, MatrixXd& Y_final, std::vector<MatrixXd>& Z_final)
{
    int n = A1.rows();
    A1 = A1.block(0, 0, n, n);
    C2 = C2.block(0, 0, n, n);
    int Z_index = 0;
    std::vector<MatrixXd> Z_g;
    MatrixXd MeanC1, MeanB1;
    batchSolveXY(C1, B1, opt, nstd1, nstd2, Z_g, MeanC1, MeanB1);

    for (int i = 0; i < Z_g.size(); i++) {
        if (Z_g[i].determinant() > 0) {
            Z_index++;
            Z_final.push_back(Z_g[i]);
        }
    }

    X_final = A1.inverse() * (B1 - C2 * Z_final[0] * B2);
    Y_final = Z_final[0] * B2;
}


/*int main() {
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
}*/

