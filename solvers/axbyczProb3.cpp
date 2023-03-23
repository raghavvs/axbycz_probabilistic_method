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
            MbMat_2(C2_m[j], Zupdate, B2inv_m[j],
                    SE3inv(Yupdate), A2_m[j], Xupdate,
                    SigB2[j], SigA2[j], B2_m[j],
                    MM[j+Ni], bb[j+Ni]);
        }

        // Solve for X and Y
        Eigen::MatrixXd M(12*(Ni+Nj), 12);
        Eigen::VectorXd b(12*(Ni+Nj));

        for (int i=0; i<Ni+Nj; ++i) {
            M.block(i*12, 0, 12, 12) = MM[i];
            b.segment(i*12, 12) = bb[i];
        }

        xi = (M.transpose() * M).ldlt().solve(M.transpose() * b);

        Eigen::Matrix4d dX;
        dX << xi(0), xi(1), xi(2), xi(3),
                xi(4), xi(5), xi(6), xi(7),
                xi(8), xi(9), xi(10), xi(11),
                0, 0, 0, 1;

        Eigen::Matrix4d dY;
        dY << xi(12), xi(13), xi(14), xi(15),
                xi(16), xi(17), xi(18), xi(19),
                xi(20), xi(21), xi(22), xi(23),
                0, 0, 0, 1;

        Xupdate = Xupdate * dX;
        Yupdate = Yupdate * dY;

        diff = metric(A1,B1,C1,Xupdate,Yupdate,Zupdate) +
               metric(A2,B2,C2,Xupdate,Yupdate,Zupdate);

        ++num;

        // Concatenate M and b matrices together
        Eigen::MatrixXd M(12*(Ni+Nj), 18);
        Eigen::VectorXd b(12*(Ni+Nj));

        for (int i=0; i<Ni+Nj; ++i) {
            M.block(i*12, 0, 12, 18) = MM[i];
            b.segment(i*12, 12) = bb[i];
        }

        // Inversion to get xi_X, xi_Y, xi_Z
        Eigen::VectorXd xi = (M.transpose() * M).ldlt().solve(M.transpose() * b);

        for (int i=0; i<Ni; ++i) {
            double diff1 = (A1_m[i] * Xupdate * B1_m[i] - Yupdate * C1_m[i] * Zupdate).norm();
        }

        for (int j=0; j<Nj; ++j) {
            double diff2 = (A2_m[j] * Xupdate * B2_m[j] - Yupdate * C2_m[j] * Zupdate).norm();
        }

        Eigen::Vector3d w_X = xi.segment(0,3);
        Eigen::Vector3d v_X = xi.segment(3,3);
        Eigen::Vector3d w_Y = xi.segment(6,3);
        Eigen::Vector3d v_Y = xi.segment(9,3);
        Eigen::Vector3d w_Z = xi.segment(12,3);
        Eigen::Vector3d v_Z = xi.segment(15,3);

        Eigen::Matrix4d X_hat;
        X_hat << skew(w_X), v_X,
                Eigen::RowVector4d::Zero();

        Eigen::Matrix4d Y_hat;
        Y_hat << skew(w_Y), v_Y,
                Eigen::RowVector4d::Zero();

        Eigen::Matrix4d Z_hat;
        Z_hat << skew(w_Z), v_Z,
                Eigen::RowVector4d::Zero();

        X_cal = Xupdate * expm(X_hat);
        Y_cal = Yupdate * expm(Y_hat);
        Z_cal = Zupdate * expm(Z_hat);

        // Update
        Xupdate = X_cal;
        Yupdate = Y_cal;
        Zupdate = Z_cal;

        ++num;

        // Error
        diff = metric(A1,B1,C1,Xupdate,Yupdate,Zupdate) +
               metric(A2,B2,C2,Xupdate,Yupdate,Zupdate);

    }

    void MbMat_1(const Ref<const MatrixXd>& A,
                 const Ref<const MatrixXd>& X,
                 const Ref<const MatrixXd>& B,
                 const Ref<const MatrixXd>& Y,
                 const Ref<const MatrixXd>& C,
                 const Ref<const MatrixXd>& Z,
                 const Ref<const MatrixXd>& SigB,
                 const Ref<const MatrixXd>& SigC,
                 Ref<MatrixXd> M,
                 Ref<VectorXd> b) {
        // Construction M and b matrices
        Eigen::Vector3d e1(1,0,0), e2(0,1,0), e3(0,0,1);

        Eigen::MatrixXd M1 = Eigen::MatrixXd::Zero(12,18);
        Eigen::VectorXd b1 = Eigen::VectorXd::Zero(12);

        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                Eigen::Vector3d Eij = A.block(0,i*4+3,3,1).cross(B.block(j*4+3,0,1,4).transpose()) +
                                      A.block(0,i*4,j+1,1).cross(B.block(j*4,j+1,i+1,3));

                M1.block(i*3+j*9  , 6+j*6  , 3 , 3) = -SigB(i,j)*skew(Eij);
                M1.block(i*3+j*9+6 , 6+j*6+3 , 3 , 3) = -SigB(i,j)*Eij;

                b1.segment(i*9+j*27   , 27) += SigB(i,j)*(C.col(j).head<3>().cross(Y.col(j).head<3>()) -
                                                          skew(C.col(j).head<3>()) * Y.col(j).head<3>() +
                                                          skew(Y.col(j).head<3>()) * C.col(j).head<>
                Z.row(i).tail<4>().transpose().cross(X.row(i)) -
                skew(Z.row(i).tail<4>()) * X.row(i) +
                skew(X.row(i)) * Z.row(i).tail<4>());

                M2.block((i-2)*9+(j-2)*27   , j-2   , 9 , 1) = SigC((i-2),(j-2))*skew(Eij);
                M2.block((i-2)*9+(j-2)*27   , j+7   , 9 , 1) = SigC((i-2),(j-2))*Eij;

                b2.segment((i-2)*27+(j-2)*81   ,81) += SigC((i-2),(j-2))*(A.col(j).tail<4>().cross(B.col(j)) -
                                                                          skew(A.col(j).tail<4>()) * B.col(j) +
                                                                          skew(B.col(j)) * A.col(j).tail<4>()) +
                                                       skew(B.col(j)) * A.col(j).tail<4>());
            }
        }

        M << M1, M2;
        b << b1, b2;
    }


    