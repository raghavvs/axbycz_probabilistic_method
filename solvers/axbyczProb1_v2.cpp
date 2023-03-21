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
*/

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "rotError.h"
#include "tranError.h"
#include "meanCov.h"
#include "so3Vec.h"

void batchSolveXY(const std::vector<Eigen::Matrix4d> &A,
                  const std::vector<Eigen::Matrix4d> &B,
                  int len,
                  bool opt,
                  double nstd_A,
                  double nstd_B,
                  std::vector<Eigen::Matrix4d> &X,
                  std::vector<Eigen::Matrix4d> &Y,
                  Eigen::MatrixXd& MeanA,
                  Eigen::MatrixXd& MeanB,
                  Eigen::MatrixXd& SigA,
                  Eigen::MatrixXd& SigB) {

    std::vector<Eigen::Matrix4d> X_candidate(8), Y_candidate(8);

    // Calculate mean and covariance for A and B
    meanCov(A, len, MeanA, SigA);
    meanCov(B, len, MeanB, SigB);

    // update SigA and SigB if nstd_A and nstd_B are known
    if (opt) {
        SigA -= nstd_A * Eigen::MatrixXd::Identity(6, 6);
        SigB -= nstd_B * Eigen::MatrixXd::Identity(6, 6);
    }

    // Calculate eigenvectors of top left 3x3 sub-matrices of SigA and SigB
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_solver_A(SigA.topLeftCorner<3, 3>());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_solver_B(SigB.topLeftCorner<3, 3>());

    auto const& VA = eig_solver_A.eigenvectors();
    auto const& VB = eig_solver_B.eigenvectors();

    // Define Q matrices
    Eigen::MatrixXd Q1, Q2, Q3, Q4;
    Q1 = Eigen::MatrixXd::Identity(3, 3);
    Q2 = (Eigen::MatrixXd(3, 3) << -1, 0, 0, 0, -1, 0, 0, 0, 1).finished();
    Q3 = (Eigen::MatrixXd(3, 3) << -1, 0, 0, 0, 1, 0, 0, 0, -1).finished();
    Q4 = (Eigen::MatrixXd(3, 3) << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished();

    Eigen::Matrix3d Rx_solved[8];

    // There are eight possibilities for Rx
    Rx_solved[0] = VA * Q1 * VB.transpose();
    Rx_solved[1] = VA * Q2 * VB.transpose();
    Rx_solved[2] = VA * Q3 * VB.transpose();
    Rx_solved[3] = VA * Q4 * VB.transpose();
    Rx_solved[4] = VA * (-Q1) * VB.transpose();
    Rx_solved[5] = VA * (-Q2) * VB.transpose();
    Rx_solved[6] = VA * (-Q3) * VB.transpose();
    Rx_solved[7] = VA * (-Q4) * VB.transpose();

    X.resize(8);
    Y.resize(8);

    for (int i = 0; i < 8; i++) {
        // block SigA and SigB to 3x3 sub-matrices
        Eigen::Matrix3d sigA_33 = SigA.block<3, 3>(0, 0);
        Eigen::Matrix3d sigB_33 = SigB.block<3, 3>(0, 3);

        Eigen::Matrix3d temp = (Rx_solved[i].transpose() * sigA_33 * Rx_solved[i]).inverse() *
                               (sigB_33 - Rx_solved[i].transpose() * SigA.block<3, 3>(0, 3) * Rx_solved[i]);

        Eigen::Vector3d tx_temp = so3Vec(temp.transpose());

        // Construct X and Y candidates
        X_candidate[i] << Rx_solved[i], -Rx_solved[i] * tx_temp, Eigen::Vector4d::Zero().transpose();
        Y_candidate[i] = MeanA * X_candidate[i] * MeanB.inverse();

        // Set the output X and Y
        X[i] = X_candidate[i];
        Y[i] = Y_candidate[i];
    }

    // Set the output MeanA, MeanB, SigA, and SigB
    for (int i = 0; i < 8; i++) {
        MeanA = MeanA * X[i] * MeanB.inverse();
    }
    MeanB = Eigen::Matrix4d::Identity();
    SigA = SigA.block<3, 3>(0, 0);
    SigB = SigB.block<3, 3>(0, 3);
}

void axbyczProb1(const Eigen::Matrix4d &A1,
                 const Eigen::Matrix4d &B1,
                 const Eigen::Matrix4d &C1,
                 const Eigen::Matrix4d &A2,
                 const Eigen::Matrix4d &B2,
                 const Eigen::Matrix4d &C2,
                 bool opt,
                 double nstd1,
                 double nstd2,
                 std::vector<Eigen::Matrix4d> &X_final,
                 std::vector<Eigen::Matrix4d> &Y_final,
                 std::vector<Eigen::Matrix4d> &Z_final) {

    //   A1 is constant with B1 and C1 free
    //   C2 is constant with A2 and B2 free

    auto A = A1.block(0, 0, 4, 4);
    auto C = C2.block(0, 0, 4, 4);

    //// ------ Solve for Z -------- //
    // A1 fixed, B1 and C1 free

    //// ------ using probability methods ------
    // calculate Z_g : all guesses of Z
    //int len = C1.size();
    int len = 8;
    std::cout << "len: " << len << std::endl;
    std::vector<Eigen::Matrix4d> Z_g(len);
    Eigen::MatrixXd MeanC, MeanB, SigA, SigB, SigC;

    std::vector<Eigen::Matrix4d> A1_vec(len), B1_vec(len), C1_vec(len),
                                A2_vec(len), B2_vec(len), C2_vec(len);

    for (int i = 0; i < len; ++i) {
        A1_vec[i] = A1.block(0, 0, 4, 4);
        B1_vec[i] = B1.block(0, 0, 4, 4);
        C1_vec[i] = C1.block(0, 0, 4, 4);
        A2_vec[i] = A2.block(0, 0, 4, 4);
        B2_vec[i] = B2.block(0, 0, 4, 4);
        C2_vec[i] = C2.block(0, 0, 4, 4);
    }

    batchSolveXY(C1_vec, B1_vec, len, opt,nstd1,nstd2,Z_g,Y_final,
                 MeanC,MeanB,SigC,SigB);

    // Keep the candidates of Z that are SE3
    // Normally there will be four Z \in SE3
    int Z_index = 0;
    std::cout << "Z_g: " << Z_g[0] << std::endl;
    std::cout << "Z_g.size(): " << Z_g.size() << std::endl;
    std::cout << "Z_g[0].size(): " << Z_g[0].size() << std::endl;
    for (int i = 0; i < Z_g.size(); ++i) {
        if (Z_g[i].determinant() > 0) {
            Z_final.push_back(Z_g[i]);
            ++Z_index;
        }
    }

    int s_Z = Z_final.size();
    std::cout << "s_Z: " << s_Z << std::endl;
    std::cout << "Z_index: " << Z_index << std::endl;
    std::cout << "Z_g: " << std::endl;
    for (int i = 0; i < Z_g.size(); ++i){
        std::cout << Z_g[i] << std::endl;
        std::cout << i << std::endl;
    }
    std::cout << "Z_final: " << std::endl;
    for (int i = 0; i < Z_final.size(); ++i){
        std::cout << Z_final[i] << std::endl;
        std::cout << i << std::endl;
    }

    std::cout << "works till here - solve for Z? - YES" << std::endl;

}

int main() {
    Eigen::Matrix4d A1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d B1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d C1 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d A2 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d B2 = Eigen::Matrix4d::Random();
    Eigen::Matrix4d C2 = Eigen::Matrix4d::Random();

    bool opt = true;
    double nstd1 = 0.1;
    double nstd2 = 0.1;

    //Eigen::Matrix4d X_final;
    std::vector<Eigen::Matrix4d> X_final;
    std::vector<Eigen::Matrix4d> Y_final;
    std::vector<Eigen::Matrix4d> Z_final;

    axbyczProb1(A1,B1,C1,A2,B2,C2,opt,nstd1,nstd2,X_final,Y_final,Z_final);

    std::cout << "Build successful?" <<std::endl;
    //std::cout << "Z_final: " << std::endl << Z_final[0] << std::endl;
    //std::cout << "X_final: " << std::endl << X_final[0] << std::endl;

    return 0;
}

/*
Output:

 Needs fixing - Z_g

s_Z: 0
Z_index: 0
Z_g:            1            0            0  5.44603e-32
0     0.835801     0.549032 -6.56859e-47
0     0.549032    -0.835801 -1.82759e-31
0            0            0            0
works till here - solve for Z? - YES
        s_X: 0
X_index: 0
works till here - solve for Z and X? - YES
        s_Y: 0
Y_index: 0
*/
