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
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include "metric.h"
#include "meanCov.h"

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
             Eigen::MatrixXd &b,
             const Eigen::MatrixXd &A,
             const Eigen::MatrixXd &X,
             const Eigen::MatrixXd &B,
             const Eigen::MatrixXd &Y,
             const Eigen::MatrixXd &C,
             const Eigen::MatrixXd &Z,
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

    M.resize(12, 18);
    M << M11, Eigen::Matrix3d::Zero(3, 3), M13, Eigen::Matrix3d::Zero(3, 3), M15, Eigen::Matrix3d::Zero(3, 3),
            M21, Eigen::Matrix3d::Zero(3, 3), M23, Eigen::Matrix3d::Zero(3, 3), M25, Eigen::Matrix3d::Zero(3, 3),
            M31, Eigen::Matrix3d::Zero(3,3) , M33, Eigen::Matrix3d::Zero(3,3) , M35, Eigen::Matrix3d::Zero (3,3),
            M41, M42,                          M43, M44,                        Eigen::Matrix3d::Zero (3,3), M46;

    // RHS
    Eigen::MatrixXd RHS = -A * X * B + Y * C * Z;

    b.resize(12, 1);
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

    M.conservativeResize(M.rows() + 36, M.cols() + 9);
    M.bottomRightCorner(36, 9) << Eigen::Matrix3d::Zero(),  M55,  M56,
            Eigen::Matrix3d::Zero(),  M65,  M66,
            Eigen::Matrix3d::Zero(),  M75,  M76,
            Eigen::Matrix3d::Zero(),  M85,  M86,
            Eigen::Matrix3d::Zero(),  M95,  M96,
            Eigen::Matrix3d::Zero(), M105, M106,
            Eigen::Matrix3d::Zero(), M115, M116,
            Eigen::Matrix3d::Zero(), M125, M126,
            Eigen::Matrix3d::Zero(), M135, M136,
            Eigen::Matrix3d::Zero(), M145, M146,
            Eigen::Matrix3d::Zero(), M155, M156,
            Eigen::Matrix3d::Zero(), M165 ,M166;

    Eigen::MatrixXd RHS2 = SE3Adinv(Z) * SigC * SE3Adinv(Z).transpose() - SigB;
    RHS2.resize(3, 12);

    b.resize(b.rows() + RHS2.size(), 1);
    b.bottomRows(RHS2.size()) = Eigen::Map<Eigen::VectorXd>(RHS2.data(), RHS2.size());
}

void MbMat_2(Eigen::MatrixXd &M,
             Eigen::MatrixXd &b,
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

    // Translation Part
    Eigen::Matrix3d M42 = -Yinv.block<3,3>(0,0) * A.block<3,3>(0,0) * X.block<3,3>(0,0);
    Eigen::Matrix3d M43 = -skew(Yinv.block<3,3>(0,0) * A.block<3,3>(0,0) * X.block<3,1>(0,3) + Yinv.block<3,3>(0,0) * A.block<3,1>(0,3) + Yinv.block<3,1>(0,3));
    Eigen::Matrix3d M44 = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d M45 = -C.block<3,3>(0,0) * Z.block<3,3>(0,0) * skew(Binv.block<3,1>(0,3));
    Eigen::Matrix3d M46 = C.block<3,3>(0,0) * Z.block<3,3>(0,0);

    M.resize(12, 18);
    M << M11, Eigen::Matrix3d::Zero(), M13, Eigen::Matrix3d::Zero(), M15, Eigen::Matrix3d::Zero(),
            M21, Eigen::Matrix3d::Zero(), M23, Eigen::Matrix3d::Zero(), M25, Eigen::Matrix3d::Zero(),
            M31, Eigen::Matrix3d::Zero(), M33, Eigen::Matrix3d::Zero(), M35, Eigen::Matrix3d::Zero(),
            Eigen::Matrix3d::Zero(), M42, M43, M44, M45, M46;

    // RHS
    Eigen::MatrixXd RHS = - C * Z * Binv + Yinv * A * X;

    b.resize(12, 1);
    b << RHS.block<3, 1>(0, 0), RHS.block<3, 1>(0, 1),
            RHS.block<3, 1>(0, 2), RHS.block<3, 1>(0, 3);

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

    M.conservativeResize(M.rows() + 36, M.cols() + 9);
    M.bottomRightCorner(36, 9) << Eigen::Matrix3d::Zero(),  M51,  M52,
            Eigen::Matrix3d::Zero(),  M61,  M62,
            Eigen::Matrix3d::Zero(),  M71,  M72,
            Eigen::Matrix3d::Zero(),  M81,  M82,
            Eigen::Matrix3d::Zero(),  M91,  M92,
            Eigen::Matrix3d::Zero(), M101, M102,
            Eigen::Matrix3d::Zero(), M111, M112,
            Eigen::Matrix3d::Zero(), M121, M122,
            Eigen::Matrix3d::Zero(), M131, M132,
            Eigen::Matrix3d::Zero(), M141, M142,
            Eigen::Matrix3d::Zero(), M151, M152,
            Eigen::Matrix3d::Zero(), M161 ,M162;

    Eigen::MatrixXd RHS2 = SE3Adinv(X) * SigA * SE3Adinv(X).transpose() - SigBinv;
    RHS2.resize(3, 12);

    b.resize(b.rows() + RHS2.size(), 1);
    b.bottomRows(RHS2.size()) = Eigen::Map<Eigen::VectorXd>(RHS2.data(), RHS2.size());
}

void axbyczProb3(const std::vector<Eigen::Matrix4d> &A1,
                 const std::vector<Eigen::Matrix4d> &B1,
                 const std::vector<Eigen::Matrix4d> &C1,
                 const std::vector<Eigen::Matrix4d> &A2,
                 const std::vector<Eigen::Matrix4d> &B2,
                 const std::vector<Eigen::Matrix4d> &C2,
                 const Eigen::Matrix4d &Xinit,
                 const Eigen::Matrix4d &Yinit,
                 const Eigen::Matrix4d &Zinit,
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
    std::vector<Eigen::Matrix4d> A1_m(Ni), B1_m(Ni), C1_m(Ni);
    std::vector<Eigen::Matrix<double, 6, 6>> SigA1(Ni), SigB1(Ni), SigC1(Ni);

    for (int i = 0; i < Ni; ++i) {
        meanCov(A1, A1_m[i], SigA1[i]);
        meanCov(B1, B1_m[i], SigB1[i]);
        meanCov(C1, C1_m[i], SigC1[i]);
    }

    // invert B2
    std::vector<Eigen::Matrix4d> B2inv(B2.size());
    for (int i = 0; i < B2.size(); ++i) {
        B2inv[i] = Eigen::MatrixXd(B2[i].rows(), B2[i].cols());
        B2inv[i] = B2[i].inverse();
    }

    std::vector<Eigen::Matrix4d> A2_m(Nj), B2_m(Nj), B2inv_m(Nj), C2_m(Nj);
    std::vector<Eigen::Matrix<double, 6, 6>> SigA2(Nj), SigB2(Nj), SigB2inv(Nj), SigC2(Nj);

    for (int j = 0; j < Nj; ++j) {
        meanCov(A2, A2_m[j], SigA2[j]);
        meanCov(B2, B2_m[j], SigB2[j]);
        meanCov(B2inv, B2inv_m[j], SigB2inv[j]);
        meanCov(C2, C2_m[j], SigC2[j]);
    }

    // Calculate M and b matrices when fixing A and C separately
    double diff = metric(A1,B1,C1,Xupdate,Yupdate,Zupdate) +
                  metric(A2,B2,C2,Xupdate,Yupdate,Zupdate);
    diff = 1;

    std::vector<Eigen::MatrixXd> MM(Ni+Nj);
    std::vector<Eigen::MatrixXd> bb(Ni+Nj);
    while (xi.norm() >= tol && diff >= tol && num <= max_num) {
        for (int i = 0; i < Ni; i++) {
            MbMat_1(MM[i], bb[i], A1_m[i], Xupdate,
                    B1_m[i], Yupdate, C1_m[i], Zupdate,
                    SigB1[i], SigC1[i]);
        }

        for (int j = 0; j < Nj; j++) {
            MbMat_2(MM[j + Ni], bb[j + Ni],
                    C2_m[j], Zupdate, B2inv_m[j],
                    SE3inv(Yupdate), A2_m[j], Xupdate,
                    SigB2[j], SigA2[j], B2_m[j]);
        }

        Eigen::MatrixXd M;
        Eigen::MatrixXd b;
        Eigen::MatrixXd M1;
        Eigen::MatrixXd M2;
        Eigen::MatrixXd M3;
        Eigen::MatrixXd M4;

        for (int k = 0; k < Ni + Nj; k++) {
            M.conservativeResize(M.rows() + MM[k].rows(), MM[k].cols());
            M.bottomRows(MM[k].rows()) = MM[k];
            b.conservativeResize(b.rows() + bb[k].rows(), b.cols() + bb[k].cols());
            b.bottomRightCorner(bb[k].rows(), bb[k].cols()) = bb[k];
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

        // Inversion to get xi_X, xi_Y, xi_Z
        //xi = (M.transpose() * M).ldlt().solve(M.transpose() * b);
        Eigen::MatrixXd xi_new = (M.transpose() * M).ldlt().solve(M.transpose() * b);

        double diff1 = 0;
        double diff2 = 0;

        for (int i = 0; i < Ni; i++) {
            diff1 = (A1_m[i] * Xupdate * B1_m[i] - Yupdate * C1_m[i] * Zupdate).norm();
        }

        for (int i = 0; i < Nj; i++) {
            diff2 = (A2_m[i] * Xupdate * B2_m[i] - Yupdate * C2_m[i] * Zupdate).norm();
        }

        Eigen::VectorXd w_X = xi_new.block(0, 0, 3, 1);
        Eigen::VectorXd v_X = xi_new.block(3, 0, 3, 1);
        Eigen::VectorXd w_Y = xi_new.block(6, 0, 3, 1);
        Eigen::VectorXd v_Y = xi_new.block(9, 0, 3, 1);
        Eigen::VectorXd w_Z = xi_new.block(12, 0, 3, 1);
        Eigen::VectorXd v_Z = xi_new.block(15, 0, 3, 1);

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

int main()
{
    srand(12345);

    std::vector<Eigen::Matrix4d> A1, B1, C1, A2, B2, C2;

    for (int i = 0; i < 10; ++i) {
        A1.emplace_back(Eigen::Matrix4d::Random());
        B1.emplace_back(Eigen::Matrix4d::Random());
        C1.emplace_back(Eigen::Matrix4d::Random());
        A2.emplace_back(Eigen::Matrix4d::Random());
        B2.emplace_back(Eigen::Matrix4d::Random());
        C2.emplace_back(Eigen::Matrix4d::Random());

        A1[i].row(3) << 0, 0, 0, 1;
        B1[i].row(3) << 0, 0, 0, 1;
        C1[i].row(3) << 0, 0, 0, 1;
        A2[i].row(3) << 0, 0, 0, 1;
        B2[i].row(3) << 0, 0, 0, 1;
        C2[i].row(3) << 0, 0, 0, 1;
    }

    Eigen::Matrix4d Xinit = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d Yinit = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d Zinit = Eigen::Matrix4d::Identity();

    Eigen::Matrix4d X_cal;
    Eigen::Matrix4d Y_cal;
    Eigen::Matrix4d Z_cal;

    int num = 5;

    axbyczProb3(A1, B1, C1, A2, B2, C2, Xinit, Yinit, Zinit, X_cal, Y_cal, Z_cal, num);

    std::cout << "Build successful? - YES" <<std::endl;

    std::cout << "X_cal: " << std::endl << X_cal << std::endl;
    std::cout << "Y_cal: " << std::endl << Y_cal << std::endl;
    std::cout << "Z_cal: " << std::endl << Z_cal << std::endl;

    return 0;
}

/*
 * Output:
 * Build successful? - YES
X_cal:
-0.533143  0.494036  0.686795  -18.2323
 0.445966  0.853957 -0.268089  0.622805
-0.718939  0.163358 -0.675605  -6.80599
        0         0         0         1
Y_cal:
  0.940418  -0.313242  -0.132265  -0.283785
-0.0757981   0.186075  -0.979608   -24.1454
  0.331465   0.931266   0.151245    1.82138
         0          0          0          1
Z_cal:
   0.992837   0.0178221    0.118138    -2.77582
 -0.0168863    0.999818 -0.00891735    -24.0595
  -0.118275  0.00685857    0.992957    -43.2199
          0           0           0           1

Build successful? - YES
X_cal:
-0.205649   0.86458  0.458486  -8.17887
-0.948059 -0.292187  0.125744  -40.0269
 0.242679 -0.408813  0.879761   96.1951
        0         0         0         1
Y_cal:
  0.607344  -0.750147  -0.261558   -46.5999
  -0.19978    0.17444  -0.964188    4.80309
  0.768909   0.637848 -0.0439189    34.7721
         0          0          0          1
Z_cal:
-0.777332  -0.59986  0.189536  -215.174
-0.624572  0.699841 -0.346601    61.419
0.0752673 -0.387803 -0.918664  -11.5001
        0         0         0         1
 */