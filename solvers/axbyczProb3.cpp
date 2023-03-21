/*
DESCRIPTION:

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

#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std; using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;

void axbyczProb3(vector<vector<Matrix4d>> &A1, vector<vector<Matrix4d>> &B1, vector<vector<Matrix4d>> &C1,
                 vector<vector<Matrix4d>> &A2, vector<vector<Matrix4d>> &B2, vector<vector<Matrix4d>> &C2,
                 Matrix4d &Xinit, Matrix4d &Yinit, Matrix4d &Zinit, Matrix4d &X_cal, Matrix4d &Y_cal,
                 Matrix4d &Z_cal) {
    // Initiation int Ni = A1.size(); int Nj = C2.size(); X_cal = Xinit; Y_cal = Yinit; Z_cal = Zinit;
    // Matrix4d Xupdate = Xinit; Matrix4d Yupdate = Yinit; Matrix4d Zupdate = Zinit; VectorXd xi(18);
    // xi.setOnes(); int max_num = 500; double tol = 1e-5; int num = 0;

    // Calculate mean and covariance of varying data
    vector<Matrix6d> SigA1(Ni);
    vector<Matrix6d> SigB1(Ni);
    vector<Matrix6d> SigC1(Ni);

    // invert B2 for i in size(B2)
    vector<vector<Matrix4d>> B2inv(Nj);
    for (int i=0; i<Nj; i++)
    {
        B2inv[i].resize(B2[i].size());
        for (int j=0; j<B2[i].size(); j++)
        {
            B2inv[i][j] =
                    B2[i][j].inverse();
        }
    }

    vector<Matrix6d> SigA2(Nj);
    vector<Matrix6d> SigB2(Nj);
    vector<Matrix6d> SigC2(Nj);

    // Calculate M and b matrices when fixing A and C separately
    double diff =
    metric(A1,B1,C1,Xupdate,Yupdate,Zupdate) + ... metric(A2,B2,C2,Xupdate,Yupdate,Zupdate);
    diff =0.0;
    while norm(xi) >= tol && diff >= tol && num <= max_num {
        for (int i=0; i<Ni; i++) { [MM[i], bb[i]] =
            MbMat_1(A_m[i],Xupdate,B_m[i],... Yupdate,C_m[i],Zupdate,... SigB(:,:,i),SigC(:,:,i)); }
            for (int j=0; j<Nj; j++) {
                [MM[j+Ni], bb[j+Ni]] =
            MbMat_2(C_m[j],Zupdate,Binv_m[j
            MbMat_2(C_m[j],Zupdate,Binv_m[j],... SE3inv(Yupdate),A_m[j],Xupdate,... SigB(:,:,j),SigA(:,:,j),B_m[j]); } // Concatenate M and b matrices together Eigen::MatrixXd M; Eigen::VectorXd b; Eigen::MatrixXd M1; Eigen::MatrixXd M2; Eigen::MatrixXd M3; Eigen::MatrixXd M4; for (int k=0; k<Ni+Nj; k++) { M = std::move(M).bottomRows(MM[k].rows()); b = std::move(b).bottomRows(bb[k].rows()); } for (int k=0; k<Ni; k++) { M1 = std::move(M1).bottomRows(MM[k].block(0,0,12,18)); M2 = std::move(M2).bottomRows(MM[k].block(12,0,9,18)); } for (int k=Ni; k<Ni+Nj; k++) { M3 = std::move(M3).bottomRows(MM[k].block(0,0,12,18)); M4 = std::move(M4).bottomRows(MM[k].block(12,0,9,18)); } // rank(M) // rank(M1) // rank(M2) // rank(M3) // rank(M4) // rank([M1;M2]) %% Inversion to get \xi_X,
        \xi_Y,
        \xi_Z xi =
        (M.transpose()*M).ldlt().solve(
        M.transpose()*b); // xi =
        pinv(M) * b;
        // xi =
        lsqminnorm(
        M,b); //
        disp(['cost is: ',
        num2str(norm(xi))]); for (int i=0;
        i<Ni;
        i++) { diff1 =
                       norm(A_m[i]*Xupdate*B_m[i]-Yupdate*C_m[i]*Zupdate,'fro'); } for (int j=0;
        j<Nj;
        j++) { diff2 =
                       norm(A_m[j]*Xupdate*B_m[j]-Yupdate*C_m[j]*Zupdate,'fro'); } w_X =
        xi.segment<3>(0); v_X =
        xi.segment<3>(3); w_Y =
        xi.segment<3>(6); v_Y =
        xi.segment<3>(9); w_Z =
        xi.segment<3>(12); v_Z =
        xi.segment<3>(15); X_hat <<
        skew(w_X),
        v_X,
        Eigen::
        Vector4d::
        Zero().transpose(); Y_hat <<
        skew(w_Y),
        v_Y,
        Eigen::
        Vector4d::
        Zero().transpose(); Z_hat <<
        skew(w_Z),
        v_Z,
        Eigen::
        Vector4d::
        Zero().transpose(); //
        X_cal = Xupdate*(Eigen










int main(){
    std::cout << "Build successful" << std::endl;
}