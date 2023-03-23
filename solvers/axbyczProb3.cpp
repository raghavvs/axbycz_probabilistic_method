/*
DESCRIPTION:

 This code defines a function called axbyczProb3 that provides iterative refinement for
 probabilistic solvers for the AXB=YCZ calibration problem. The function takes in the

 calibration data in the form of cell arrays A1,B1,C1 and A2,B2,C2 that store data collected
 while fixing A and C at different poses respectively. The initial guesses for X, Y, and Z
 matrices are also provided as Xinit,Yinit,Zinit. The output of the function are the calibrated
 X_cal,Y_cal,Z_cal matrices. The function uses the mean and covariance of the calibration data
 to calculate matrices M and b that are used to solve for X_cal,Y_cal,Z_cal matrices using an
 iterative solver. The function returns the number of iterations used to solve the problem.

 This function provides iterative refinement for probabilistic solvers for
  the AXB=YCZ calibration problem.
 It solves the system of equations:
  (1)  A_i X \bar{B_i} = Y \bar{B_i} Z
  (2)  \Sigma^1_{B_i} = R_Z^T \Sigma^1_{C_i} R_Z
  (3)  C_j Z \bar{B_j}^{-1} = Y^{-1} \bar{A_j} X
  (4)  R_{\bar{B_j}}^T \Sigma^1_{B_j} R_{\bar{B_j}}^T = R_X^T \Sigma^1_{A_j} R_X
 where X = X_init(I+\xi_X), Y = Y_init(I+\xi_Y), Z = Z_init(I+\xi_Z)

 Inputs:
   A1,B1,C1: Cell arrays, which stores when fixing A1 at different poses;
   A2,B2,C2: Cell arrays, which stores when fixing C2 at different poses;
   Xinit,Yinit,Zinit: Initial guesses of X,Y,Z matrices.
 Outputs:
   X_cal,Y_cal,Z_cal: Calibrated X,Y,Z matrices

 Updates:
   Used the whole covariance matrix

Author: Sipu Ruan, ruansp@jhu.edu, November 2017 (MATLAB Version)
*/

#include <Eigen/Dense>
#include <vector>

void axbyczProb3(const std::vector<Eigen::MatrixXd> &A1,
                 const std::vector<Eigen::MatrixXd> &B1,
                 const std::vector<Eigen::MatrixXd> &C1,
                 const std::vector<Eigen::MatrixXd> &A2,
                 const std::vector<Eigen::MatrixXd> &B2,
                 const std::vector<Eigen::MatrixXd> &C2,
                 Eigen::Matrix4d Xinit,
                 Eigen::Matrix4d Yinit,
                 Eigen::Matrix4d Zinit,
                 Eigen::Matrix4d &X_cal,
                 Eigen::Matrix4d &Y_cal,
                 Eigen::Matrix4d &Z_cal,
                 int& num) {
    // Initiation
    int Ni = A1.size();
    int Nj = C2.size();
    X_cal = Xinit;
    Y_cal = Yinit;
    Z_cal = Zinit;
    Eigen::Matrix4d Xupdate = Xinit;
    Eigen::Matrix4d Yupdate = Yinit;
    Eigen::Matrix4d Zupdate = Zinit;
    Eigen:VectorXd xi = Eigen:VectorXd(18);
    xi.setOnes();

    int max_num = 500;
    double tol = 1e-5;

    // Calculate mean and covariance of varying data
    std:vector<Eigen:MatrixXd> A1_m(Ni), SigA1(Ni), B1_m(Ni), SigB1(Ni), C1_m(Ni), SigC1(Ni);
    for (int i=0; i<Ni; ++i) {
        meanCov(A1[i], A1_m[i], SigA1[i]);
        meanCov(B1[i], B1_m[i], SigB1[i]);
        meanCov(C1[i], C1_m[i], SigC1[i]);
    }

    // invert B2
    std::vector<std::vector<Eigen::MatrixXd>> B2inv(B2.size());
    for (int i=0; i<B2.size(); ++i) {
        B2inv[i].resize(B2[i].size());
        for (int j=0; j<B2[i].size(); ++j) {
            B2inv[i][j] = B2[i][j].inverse();
        }
    }

    std::vector<Eigen::MatrixXd> A2_m(Nj), SigA2(Nj), B2_m(Nj), SigB2(Nj), C2_m(Nj), SigC2(Nj);
    std::vector<Eigen::MatrixXd> B2inv_m(Nj), SigB2inv(Nj);

    for (int j=0; j<Nj; ++j) {
        meanCov(A2[j], A2_m[j], SigA2[j]);
        meanCov(B2[j], B2_m[j], SigB2[j]);
        meanCov(B2inv[j], B2inv_m[j], SigB2inv[j]);
        meanCov(C2[j], C2_m[j], SigC2[j]);
    }

    // Calculate M and b matrices when fixing A and C separately
    double diff = metric(A1,B1,C1,Xupdate,Yupdate,Zupdate) +
                  metric(A2,B2,C2,Xupdate,Yupdate,Zupdate);
    diff = 1;

    while (xi.norm() >= tol && diff >= tol && num <= max_num) {
        std::vector<Eigen::MatrixXd> MM(Ni+Nj), bb(Ni+Nj);

        for (int i=0; i<Ni; ++i) {
            MbMat_1(A1_m[i], Xupdate, B1_m[i],
                    Yupdate, C1_m[i], Zupdate,
                    SigB1[i], SigC1[i],
                    MM[i], bb[i]);
        }

        for (int j=0; j<Nj; ++j) {
    }