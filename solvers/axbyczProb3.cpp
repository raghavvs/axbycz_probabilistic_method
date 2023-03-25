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

Author: Sipu Ruan, ruansp@jhu.edu, November 2017 (MATLAB Version)

 The functions SE3inv, SE3Ad, and SE3Adinv are used to perform operations on elements of the
 Special Euclidean group SE(3), which represents rigid body transformations in 3D space.

SE3inv computes the inverse of an element of SE(3). Given a transformation matrix X that
 represents a rotation and translation in 3D space, SE3inv(X) returns the transformation
 matrix that “undoes” the transformation represented by X.

SE3Ad computes the adjoint representation of an element of SE(3). Given a transformation
 matrix X, SE3Ad(X) returns a 6x6 matrix that can be used to transform spatial motion vectors
 (e.g., twists) from one coordinate frame to another.

SE3Adinv computes the inverse of the adjoint representation of an element of SE(3). Given a
 transformation matrix X, SE3Adinv(X) returns a 6x6 matrix that can be used to transform
 spatial motion vectors from one coordinate frame to another, in the opposite direction as SE3Ad(X).
*/

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "meanCov.h"
#include "metric.h"

// Define a function that takes a vector of size 3 and returns a 3x3 matrix
Eigen::Matrix3d skew(Eigen::Vector3d v) {
    Eigen::Matrix3d M;

    // Fill in the elements of the matrix according to the formula
    M << 0, -v(2), v(1),
            v(2), 0, -v(0),
            -v(1), v(0), 0;

    return M;
}

Eigen::Matrix4d SE3inv(const Eigen::Matrix4d& X) {
    Eigen::Matrix4d invX;
    invX << X.block<3,3>(0,0).transpose(), -X.block<3,3>(0,0).transpose() * X.block<3,1>(0,3),
            0, 0, 0, 1;
    return invX;
}

Eigen::Matrix<double, 6, 6> SE3Adinv(const Eigen::Matrix4d& X) {
    Eigen::Matrix3d R = X.block<3,3>(0,0);
    Eigen::Vector3d t = X.block<3,1>(0,3);

    Eigen::Matrix<double, 6, 6> A;
    A << R.transpose(), Eigen::Matrix3d::Zero(),
            -R.transpose() * Eigen::Matrix3d(Eigen::AngleAxisd(t.norm(), t.normalized())) * R, R.transpose();
    return A;
}

Eigen::Matrix<double, 6, 6> SE3Ad(const Eigen::Matrix4d& X) {
    Eigen::Matrix3d R = X.block<3,3>(0,0);
    Eigen::Vector3d t = X.block<3,1>(0,3);

    Eigen::Matrix<double, 6, 6> A;
    A << R, Eigen::Matrix3d::Zero(),
            skew(t) * R, R;
    return A;
}

void MbMat_1(Eigen::MatrixXd &M,
             Eigen::VectorXd &b,
             const Eigen::Matrix4d &A,
             const Eigen::Matrix4d &X,
             const Eigen::Matrix4d &B,
             const Eigen::Matrix4d &Y,
             const Eigen::Matrix4d &C,
             const Eigen::Matrix4d &Z,
             const Eigen::MatrixXd &SigB,
             const Eigen::MatrixXd &SigC) {
    // Construction M and b matrices
    Eigen::Vector3d e1(1, 0, 0);
    Eigen::Vector3d e2(0, 1, 0);
    Eigen::Vector3d e3(0, 0, 1);

    // AXB = YCZ
    // Rotation part
    Eigen::Matrix3d M11 = -A.block<3,3>(0,0) * X.block<3,3>(0,0) * skew(B.block<3,3>(0,0)*e1);
    Eigen::Matrix3d M13 = Y.block<3,3>(0,0) * skew(C.block<3,3>(0,0)*Z.block<3,3>(0,0)*e1);
    Eigen::Matrix3d M15 = Y.block<3,3>(0,0) * C.block<3,3>(0,0) * Z.block<3,3>(0,0) * skew(e1);

    Eigen::Matrix3d M21 = -A.block<3,3>(0,0) * X.block<3,3>(0,0) * skew(B.block<3,3>(0,0)*e2);
    Eigen::Matrix3d M23 = Y.block<3,3>(0,0) * skew(C.block<3,3>(0,0)*Z.block<3,3>(0,0)*e2);
    Eigen::Matrix3d M25 = Y.block<3,3>(0,0) * C.block<3,3>(0,0) * Z.block<3,3>(0,0) * skew(e2);

    Eigen::Matrix3d M31 = -A.block<3,3>(0,0) * X.block<3,3>(0,0) * skew(B.block<3,3>(0,0)*e3);
    Eigen::Matrix3d M33 = Y.block<3,3>(0,0) * skew(C.block<3,3>(0,0)*Z.block<3,3>(0,0)*e3);
    Eigen::Matrix3d M35 = Y.block<3,3>(0,0) * C.block<3,3>(0,0) * Z.block<3,3>(0,0) * skew(e3);

    // Translation part
    Eigen::Matrix3d M41 = -A.block<3,3>(0,0) * X.block<3,3>(0,0) * skew(B.block<3,1>(0,3));
    Eigen::Matrix3d M42 = A.block<3,3>(0,0) * X.block<3,3>(0,0);
    Eigen::Matrix3d M43 = Y.block<3, 3>(0, 0) * skew(C.block<3, 3>(0, 0) * Z.block<3,1>(0,3) + C.block<3,1>(0,3));
    Eigen::Matrix3d M44 = -Y.block<3,3>(0,0);
    Eigen::Matrix3d M46 = -Y.block<3,3>(0,0) * C.block<3,3>(0,0) * Z.block<3,3>(0,0);

    M.resize(12, 12);
    M << M11, Eigen::Matrix3d::Zero(3, 3), M13, Eigen::Matrix3d::Zero(3, 3), M15, Eigen::Matrix3d::Zero(3, 3),
            M21, Eigen::Matrix3d::Zero(3, 3), M23, Eigen::Matrix3d::Zero(3, 3), M25, Eigen::Matrix3d::Zero(3, 3),
            M31, Eigen::Matrix3d::Zero(3,3) , M33, Eigen::Matrix3d::Zero(3,3) , M35, Eigen::Matrix3d::Zero (3,3),
            M41, M42,                          M43, M44,                        Eigen::Matrix3d::Zero (3,3), M46;

    // RHS
    Eigen::MatrixXd RHS = -A * X * B + Y * C * Z;

    b << RHS.block<3, 1>(0, 0), RHS.block<3, 1>(0, 1), RHS.block<3, 1>(0, 2), RHS.block<3, 1>(0, 3);

    // SigBi = Ad^{-1}(Z) * SigCi * Ad^{-T}(Z)
    // First block
    Eigen::Matrix3d M55 = -skew(SigB.block<3,1>(0,0)) + SigB.block<3,3>(0,0) * skew(e1);
    Eigen::Matrix3d M56 = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d M65 = -skew(SigB.block<3,1>(0,1)) + SigB.block<3,3>(0,0) * skew(e2);
    Eigen::Matrix3d M66 = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d M75 = -skew(SigB.block<3,1>(0,2)) + SigB.block<3,3>(0,0) * skew(e3);
    Eigen::Matrix3d M76 = Eigen::Matrix3d::Zero();

    // Second block
    Eigen::Matrix3d M85 = -skew(SigB.block<3,1>(0,3)) + SigB.block<3,3>(0,3) * skew(e1);
    Eigen::Matrix3d M86 = SigB.block<3,3>(0,0) * skew(e1);
    Eigen::Matrix3d M95 = -skew(SigB.block<3,1>(0,4)) + SigB.block<3,3>(0,3) * skew(e2);
    Eigen::Matrix3d M96 = SigB.block<3,3>(0,0) * skew(e2);
    Eigen::Matrix3d M105 = -skew(SigB.block<3,1>(0,5)) + SigB.block<3,3>(0,3) * skew(e3);
    Eigen::Matrix3d M106 = SigB.block<3,3>(0,0) * skew(e3);

    // Third block
    Eigen::Matrix3d M115 = -skew(SigB.block<3,1>(3,0)) + SigB.block<3,3>(3,0) * skew(e1);
    Eigen::Matrix3d M116 = -skew(SigB.block<3,1>(0,0));
    Eigen::Matrix3d M125 = -skew(SigB.block<3,1>(3,1)) + SigB.block<3,3>(3,0) * skew(e2);
    Eigen::Matrix3d M126 = -skew(SigB.block<3,1>(0,1));
    Eigen::Matrix3d M135 = -skew(SigB.block<3,1>(3,2)) + SigB.block<3,3>(3,0) * skew(e3);
    Eigen::Matrix3d M136 = -skew(SigB.block<3,1>(0,2));

    // Fourth block
    Eigen::Matrix3d M145 = -skew(SigB.block<3,1>(3,3)) + SigB.block<3,3>(3,3) * skew(e1);
    Eigen::Matrix3d M146 = -skew(SigB.block<3,1>(0,3)) + SigB.block<3,3>(3,0) * skew(e1);
    Eigen::Matrix3d M155 = -skew(SigB.block<3,1>(3,4)) + SigB.block<3,3>(3,3) * skew(e2);
    Eigen::Matrix3d M156 = -skew(SigB.block<3,1>(0,4)) + SigB.block<3,3>(3,0) * skew(e2);
    Eigen::Matrix3d M165 = -skew(SigB.block<3,1>(3,5)) + SigB.block<3,3>(3,3) * skew(e3);
    Eigen::Matrix3d M166 = -skew(SigB.block<3,1>(0,5)) + SigB.block<3,3>(3,0) * skew(e3);

    M.conservativeResize(M.rows() + 36, Eigen::NoChange);
    M.block(M.rows() - 36, 0, 3, 12) << M55, M56;
    M.block(M.rows() - 33, 0, 3, 12) << M65, M66;
    M.block(M.rows() - 30, 0, 3, 12) << M75, M76;
    M.block(M.rows() - 27, 0, 3, 12) << M85, M86;
    M.block(M.rows() - 24, 0, 3, 12) << M95, M96;
    M.block(M.rows() - 21, 0, 3, 12) << M105,M106;
    M.block(M.rows() - 18, 0, 3, 12) << M115,M116;
    M.block(M.rows() - 15, 0, 3, 12) << M125,M126;
    M.block(M.rows() - 12, 0, 3, 12) << M135,M136;
    M.block(M.rows() - 9 , 0, 3, 12) << M145,M146;
    M.block(M.rows() - 6 , 0, 3, 12) << M155,M156;

    Eigen::MatrixXd RHS2 = SE3Adinv(Z) * SigC * SE3Adinv(Z).transpose() - SigB;
    RHS2.resize(3, 12);

    b.conservativeResize(b.size() + 36);
    b.tail(36) = Eigen::Map<Eigen::VectorXd>(RHS2.data(), 36);
}

void MbMat_2(Eigen::MatrixXd &M,
             Eigen::VectorXd &b,
             const Eigen::MatrixXd &C,
             const Eigen::MatrixXd &Z,
             const Eigen::MatrixXd &Binv,
             const Eigen::MatrixXd &Yinv,
             const Eigen::MatrixXd &A,
             const Eigen::MatrixXd &X,
             const Eigen::MatrixXd &SigB,
             const Eigen::MatrixXd &SigA,
             const Eigen::MatrixXd &B){
    // Construction M and b matrices
    Eigen::Vector3d e1(1,0,0), e2(0,1,0), e3(0,0,1);
    Eigen::Matrix3d Binv3 = B.topLeftCorner<3,3>().inverse();
    Eigen::MatrixXd SigBinv = SE3Ad(B) * SigB * SE3Ad(B).transpose();

    // CZB^{-1} = Y^{-1}AX
    // Rotation part
    Eigen::Matrix3d M11 = Yinv.topLeftCorner<3,3>() * A.topLeftCorner<3,3>() * X.topLeftCorner<3,3>() * skew(e1);
    Eigen::Matrix3d M13 = -skew(Yinv.topLeftCorner<3,3>() * A.topLeftCorner<3,3>() * X.topLeftCorner<3,3>() * e1);
    Eigen::Matrix3d M15 = -C.topLeftCorner<3,3>() * Z.topLeftCorner<3,3>() * skew(Binv3*e1);

    Eigen::Matrix3d M21 = Yinv.topLeftCorner<3,3>() * A.topLeftCorner<3,3>() * X.topLeftCorner<3,3>() * skew(e2);
    Eigen::Matrix3d M23 = -skew(Yinv.topLeftCorner<3,3>() * A.topLeftCorner<3,3>() * X.topLeftCorner<3,3>() * e2);
    Eigen::Matrix3d M25 = -C.topLeftCorner<3,3>() * Z.topLeftCorner<3,3>() * skew(Binv3*e2);

    Eigen::Matrix3d M31 = Yinv.topLeftCorner<3,3>() * A.topLeftCorner<3,3>() * X.topLeftCorner<3,3>() * skew(e3);
    Eigen::Matrix3d M33 = -skew(Yinv.topLeftCorner<3,3>() * A.topLeftCorner<3,3>() * X.topLeftCorner<3,3>() * e3);
    Eigen::Matrix3d M35 = -C.topLeftCorner<3,3>() * Z.topLeftCorner<3,3>() * skew(Binv3*e3);

    // Translation part
    Eigen::Matrix3d M42 = -Yinv.topLeftCorner<3,3>() * A.topLeftCorner<3,3>() * X.topLeftCorner<3,3>();
    Eigen::Matrix3d M43 = -skew(Yinv.topLeftCorner<3,3>()*A.topLeftCorner<3,3>()*X.col(3) + Yinv.topLeftCorner<3,3>()*A.col(3) + Yinv.col(3));
    Eigen::Matrix3d M44 = Eigen::MatrixXd::Identity(3, 3);
    Eigen::Matrix3d M45 = -C.topLeftCorner<3,3>() * Z.topLeftCorner<3,3>() * skew(Binv.col(3));

    Eigen::Matrix3d M46 = C.topLeftCorner<3,3>() * Z.topLeftCorner<3,3>();

    M.resize(12, 12);
    M << M11, Eigen::Matrix3d::Zero(), M13, Eigen::Matrix3d::Zero(), M15, Eigen::Matrix3d::Zero(),
            M21, Eigen::Matrix3d::Zero(), M23, Eigen::Matrix3d::Zero(), M25, Eigen::Matrix3d::Zero(),
            M31, Eigen::Matrix3d::Zero(), M33, Eigen::Matrix3d::Zero(), M35, Eigen::Matrix3d::Zero(),
            Eigen::Matrix3d::Zero(), M42, M43, M44, M45, M46;

    // RHS
    Eigen::MatrixXd RHS = - C * Z * Binv + Yinv * A * X;

    b.resize(12);
    b << RHS.block<3,1>(0,0), RHS.block<3,1>(0,1), RHS.block<3,1>(0,2), RHS.block<3,1>(0,3);

    // SigBi^{-1} = Ad^{-1}(X) * SigAi * Ad^{-T}(X)
    // First block
    Eigen::Matrix3d M51 = -skew(SigBinv.block<3,1>(0,0)) + SigBinv.block<3,3>(0,0) * skew(e1);
    Eigen::Matrix3d M52 = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d M61 = -skew(SigBinv.block<3,1>(0,1)) + SigBinv.block<3,3>(0,0) * skew(e2);
    Eigen::Matrix3d M62 = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d M71 = -skew(SigBinv.block<3,1>(0,2)) + SigBinv.block<3,3>(0,0) * skew(e3);
    Eigen::Matrix3d M72 = Eigen::Matrix3d::Zero();

    // Second block
    Eigen::Matrix3d M81 = -skew(SigBinv.block<3,1>(0,3)) + SigBinv.block<3,3>(0,3) * skew(e1);
    Eigen::Matrix3d M82 = SigBinv.block<3,3>(0,0) * skew(e1);
    Eigen::Matrix3d M91 = -skew(SigBinv.block<3,1>(0,4)) + SigBinv.block<3,3>(0,3) * skew(e2);
    Eigen::Matrix3d M92 = SigBinv.block<3,3>(0,0) * skew(e2);
    Eigen::Matrix3d M101 = -skew(SigBinv.block<3,1>(0,5)) + SigBinv.block<3,3>(0,3) * skew(e3);
    Eigen::Matrix3d M102 = SigBinv.block<3,3>(0,0) * skew(e3);

    // Third block
    Eigen::Matrix3d M111 = -skew(SigBinv.block<3,1>(3,0)) + SigBinv.block<3,3>(3,0) * skew(e1);
    Eigen::Matrix3d M112 = -skew(SigBinv.block<3,1>(0,0));
    Eigen::Matrix3d M121 = -skew(SigBinv.block<3,1>(3,1)) + SigBinv.block<3,3>(3,0) * skew(e2);
    Eigen::Matrix3d M122 = -skew(SigBinv.block<3,1>(0,1));
    Eigen::Matrix3d M131 = -skew(SigBinv.block<3,1>(3,2)) + SigBinv.block<3,3>(3,0) * skew(e3);
    Eigen::Matrix3d M132 = -skew(SigBinv.block<3,1>(0,2));

    // Fourth block
    Eigen::Matrix3d M141 = -skew(SigBinv.block<3,1>(3,3)) + SigBinv.block<3,3>(3,3) * skew(e1);
    Eigen::Matrix3d M142 = -skew(SigBinv.block<3,1>(0,3)) + SigBinv.block<3,3>(3,0) * skew(e1);
    Eigen::Matrix3d M151 = -skew(SigBinv.block<3,1>(3,4)) + SigBinv.block<3,3>(3,3) * skew(e2);
    Eigen::Matrix3d M152 = -skew(SigBinv.block<3,1>(0,4)) + SigBinv.block<3,3>(3,0) * skew(e2);
    Eigen::Matrix3d M161 = -skew(SigBinv.block<3,1>(3,5)) + SigBinv.block<3,3>(3,3) * skew(e3);
    Eigen::Matrix3d M162 = -skew(SigBinv.block<3,1>(0,5)) + SigBinv.block<3,3>(3,0) * skew(e3);

    M.conservativeResize(M.rows()+36, Eigen::NoChange);
    M.block(M.rows()-36,  0, 3, 12) << M51, M52;
    M.block(M.rows()-33,  0, 3, 12) << M61, M62;
    M.block(M.rows()-30,  0, 3, 12) << M71, M72;
    M.block(M.rows()-27,  0, 3, 12) << M81, M82;
    M.block(M.rows()-24,  0, 3, 12) << M91, M92;
    M.block(M.rows()-21,  0, 3, 12) << M101,M102;
    M.block(M.rows()-18,  0, 3, 12) << M111,M112;
    M.block(M.rows()-15,  0, 3, 12) << M121,M122;
    M.block(M.rows()-12,  0, 3, 12) << M131,M132;
    M.block(M.rows()-9 ,  0, 3, 12) << M141,M142;
    M.block(M.rows()-6 ,  0, 3, 12) << M151,M152;

    Eigen::MatrixXd RHS2 = SE3Adinv(X) * SigA * SE3Adinv(X).transpose() - SigBinv;
    RHS2.resize(3, 12);

    b = Eigen::VectorXd::Map(RHS2.data(), RHS2.size());
}

void axbyczProb3(const std::vector<Eigen::Matrix4d> &A1,
                 const std::vector<Eigen::Matrix4d> &B1,
                 const std::vector<Eigen::Matrix4d> &C1,
                 const std::vector<Eigen::Matrix4d> &A2,
                 const std::vector<Eigen::Matrix4d> &B2,
                 const std::vector<Eigen::Matrix4d> &C2,
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
    Eigen::VectorXd xi = Eigen::VectorXd(18);
    xi.setOnes();

    int max_num = 500;
    double tol = 1e-5;

    // Calculate mean and covariance of varying data
    Eigen::MatrixXd A1_m, B1_m, C1_m, SigA1, SigB1, SigC1;
    for (int i = 0; i < Ni; i++) {
        meanCov(A1, Ni, A1_m, SigA1);
        meanCov(B1, Ni, B1_m, SigB1);
        meanCov(C1, Ni, C1_m, SigC1);
    }

    // invert B2
    std::vector<Eigen::Matrix4d> B2inv(B2.size());
    for (int i = 0; i < B2.size(); ++i) {
        B2inv[i] = Eigen::MatrixXd(B2[i].rows(), B2[i].cols());
        for (int j = 0; j < B2[i].cols(); ++j) {
            B2inv[i].col(j) = B2[i].col(j).inverse();
        }
    }

    Eigen::MatrixXd A2_m, B2_m, B2inv_m, C2_m, SigA2, SigB2, SigB2inv, SigC2;
    for (int j = 0; j < Nj; j++) {
        meanCov(A2, Nj, A2_m, SigA2);
        meanCov(B2, Nj, B2_m, SigB2);
        meanCov(B2inv, Nj, B2inv_m, SigB2inv);
        meanCov(C2, Nj, C2_m, SigC2);
    }

    // Calculate M and b matrices when fixing A and C separately
    double diff = metric(A1,B1,C1,Xupdate,Yupdate,Zupdate) +
                  metric(A2,B2,C2,Xupdate,Yupdate,Zupdate);
    diff = 1;

    while (xi.norm() >= tol && diff >= tol && num <= max_num) {
        std::vector<Eigen::MatrixXd> MM(Ni+Nj);
        std::vector<Eigen::VectorXd> bb(Ni+Nj);
        for (int i = 0; i < Ni; i++) {
            MbMat_1(MM[i], bb[i], A1_m.block(0, 0, 4, 4 * i), Xupdate,
                    B1_m.block(0, 0, 4, 4 * i),
                    Yupdate, C1_m.block(0, 0, 4, 4 * i), Zupdate,
                    SigB1.block(0, 0, SigB1.rows(), SigB1.cols() * i),
                    SigC1.block(0, 0, SigC1.rows(), SigC1.cols() * i));
            MbMat_1(MM[i], bb[i], A1_m.block<4,4>(0, 4 * i), Xupdate,
                    B1_m.block<4,4>(0, 4 * i),
                    Yupdate, C1_m.block<4,4>(0, 4 * i), Zupdate,
                    SigB1.block(0, SigB1.cols() * i, SigB1.rows(), SigB1.cols()),
                    SigC1.block(0, SigC1.cols() * i, SigC1.rows(), SigC1.cols()));
        }

        for (int j = 0; j < Nj; j++) {
            MbMat_2(MM[j + Ni], bb[j + Ni],
                    C2_m.block(0, 0, C2_m.rows(), C2_m.cols() * j), Zupdate,
                    B2inv_m.block(0, 0, B2inv_m.rows(), B2inv_m.cols() * j),
                    SE3inv(Yupdate),
                    A2_m.block(0, 0, A2_m.rows(),A2_m.cols() * j), Xupdate,
                    SigB2.block(0, 0, SigB2.rows(), SigB2.cols() * j),
                    SigA2.block(0, 0, SigA2.rows(), SigA2.cols() * j),
                    B2_m.block(0, 0, B2_m.rows(), B2_m.cols() * j));
        }

        Eigen::MatrixXd M;
        Eigen::VectorXd b;
        Eigen::MatrixXd M1;
        Eigen::MatrixXd M2;
        Eigen::MatrixXd M3;
        Eigen::MatrixXd M4;

        for (int k = 0; k < Ni + Nj; k++) {
            M.conservativeResize(M.rows() + MM[k].rows(), MM[k].cols());
            M.bottomRows(MM[k].rows()) = MM[k];
            b.conservativeResize(b.size() + bb[k].size());
            b.tail(bb[k].size()) = bb[k];
        }

        for (int k = 0; k < Ni; k++) {
            M1.conservativeResize(M1.rows() + MM[k].block(0, 0, 12, MM[k].cols()).rows(), MM[k].cols());
            M1.bottomRows(MM[k].block(0, 0, 12, MM[k].cols()).rows()) = MM[k].block(0, 0, 12, MM[k].cols());
            M2.conservativeResize(M2.rows() + MM[k].block(12, 0, 9, MM[k].cols()).rows(), MM[k].cols());
            M2.bottomRows(MM[k].block(12, 0, 9, MM[k].cols()).rows()) = MM[k].block(12, 0, 9, MM[k].cols());
        }

        for (int k = Ni; k < Ni + Nj; k++) {
            M3.conservativeResize(M3.rows() + MM[k].block(0, 0, 12, MM[k].cols()).rows(), MM[k].cols());
            M3.bottomRows(MM[k].block(0, 0, 12, MM[k].cols()).rows()) = MM[k].block(0, 0, 12, MM[k].cols());
            M4.conservativeResize(M4.rows() + MM[k].block(12, 0, 9, MM[k].cols()).rows(), MM[k].cols());
            M4.bottomRows(MM[k].block(12, 0, 9, MM[k].cols()).rows()) = MM[k].block(12, 0, 9, MM[k].cols());
        }

        Eigen::VectorXd xi = (M.transpose() * M).ldlt().solve(M.transpose() * b);

        double diff1 = 0;
        double diff2 = 0;

        for (int i = 0; i < Ni; i++) {
            diff1 = (A1_m.block(0, 0, 4, 4 * i) * Xupdate * B1_m.block(0, 0, 4, 4 * i) - Yupdate * C1_m.block(0, 0, 4, 4 * i) * Zupdate).norm();
        }

        for (int i = 0; i < Nj; i++) {
            diff2 = (A2_m.block(0, 0, A2_m.rows(), A2_m.cols() * i) * Xupdate * B2_m.block(0, 0, B2_m.rows(), B2_m.cols() * i) - Yupdate * C2_m.block(0, 0, C2_m.rows(), C2_m.cols() * i) * Zupdate).norm();
        }

        Eigen::Vector3d w_X = xi.head(3);
        Eigen::Vector3d v_X = xi.segment(3, 3);
        Eigen::Vector3d w_Y = xi.segment(6, 3);
        Eigen::Vector3d v_Y = xi.segment(9, 3);
        Eigen::Vector3d w_Z = xi.segment(12, 3);
        Eigen::Vector3d v_Z = xi.tail(3);

        Eigen::Matrix4d X_hat;
        X_hat << skew(w_X), v_X,
                Eigen::RowVector4d::Zero();

        Eigen::Matrix4d Y_hat;
        Y_hat << skew(w_Y), v_Y,
                Eigen::RowVector4d::Zero();

        Eigen::Matrix4d Z_hat;
        Z_hat << skew(w_Z), v_Z,
                Eigen::RowVector4d::Zero();

        X_cal = Xupdate * (X_hat).exp();
        Y_cal = Yupdate * (Y_hat).exp();
        Z_cal = Zupdate * (Z_hat).exp();

        // Update
        Xupdate = X_cal;
        Yupdate = Y_cal;
        Zupdate = Z_cal;

        ++num;

        // Error
        diff = metric(A1,B1,C1,Xupdate,Yupdate,Zupdate) +
               metric(A2,B2,C2,Xupdate,Yupdate,Zupdate);

    }
}