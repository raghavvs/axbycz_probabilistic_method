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
#include <cstdlib>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "rotError.h"
#include "tranError.h"
#include "so3Vec.h"

void meanCov(const std::vector<Eigen::Matrix4d> &X, int N, Eigen::MatrixXd &Mean,
             Eigen::MatrixXd &Cov) {
    Mean = Eigen::Matrix4d::Identity();
    Cov = Eigen::Matrix<double, 6, 6>::Zero();

    // Initial approximation of Mean
    Eigen::Matrix4d sum_se = Eigen::Matrix4d::Zero();
    for (int i = 0; i < N; i++) {
        sum_se += X[i].log();
    }
    Mean = ((1.0 / N) * sum_se).exp();

    //sum_se is fine

    // Iterative process to calculate the true Mean
    Eigen::Matrix4d diff_se = Eigen::Matrix4d::Ones();
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num) {
        diff_se = Eigen::Matrix4d::Zero();
        for (int i = 0; i < N; i++) {
            diff_se += (Mean.inverse() * X[i]).log();
        }
        Mean *= ((1.0 / N) * diff_se).exp();
        count++;
    }
    /*std::cout << "diff_se: " << std::endl << diff_se << std::endl;
    std::cout << "X[0]: " << std::endl << X[0] << std::endl;
    std::cout << "Mean.inverse(): " << std::endl << Mean.inverse() << std::endl;
    std::cout << "Mean.inverse() * X[0]: " << std::endl << Mean.inverse() * X[0] << std::endl;*/
    //std::cout << "(Mean.inverse() * X[0]).log(): " << std::endl << (Mean.inverse() * X[0]).log() << std::endl;
    // something is causing this - (Mean.inverse() * X[0]).log() - to go berserk

    // Covariance
    for (int i = 0; i < N; i++) {
        diff_se = (Mean.inverse() * X[i]).log();
        //std::cout << "diff_se: " << std::endl << diff_se << std::endl;
        Eigen::VectorXd diff_vex(6);
        diff_vex << Eigen::Map<Eigen::Vector3d>(diff_se.block<3,3>(0,0).data()), diff_se.block<3,1>(0,3);
        Cov += diff_vex * diff_vex.transpose();
    }
    Cov /= N;

    //std::cout << "SigA: " << std::endl << Cov << std::endl;
}

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

    /*std::cout << "A: " << A[0] << std::endl;
    std::cout << "B: " << B[0] << std::endl;
    std::cout << "MeanA: " << MeanA << std::endl;
    std::cout << "MeanB: " << MeanB << std::endl;*/
    std::cout << "SigA: " << std::endl << SigA << std::endl;
    std::cout << "SigB: " << std::endl << SigB << std::endl;

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

    std::cout << "C1_vec: " << std::endl << C1_vec[0] << std::endl;
    std::cout << "B1_vec: " << std::endl << B1_vec[0] << std::endl;

    batchSolveXY(C1_vec, B1_vec, len, opt,nstd1,nstd2,Z_g,Y_final,
                 MeanC,MeanB,SigC,SigB);

    // Keep the candidates of Z that are SE3
    // Normally there will be four Z \in SE3
    double threshold = 1e-5;
    for (int i = 0; i < len; ++i){
        Z_g[i] = Z_g[i].unaryExpr([threshold](double x) { return std::abs(x) < threshold ? 0 : x; });
    }
    int Z_index = 0;
    std::cout << "Z_g: " << Z_g[0] << std::endl;
    std::cout << "Z_g.size(): " << Z_g.size() << std::endl;
    std::cout << "Z_g[0].size(): " << Z_g[0].size() << std::endl;
    std::cout << "Z_g[0].determinant(): " << Z_g[0].determinant() << std::endl;
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
    srand(12345);
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


Output new:

len: 8
C1_vec:  0.964414  0.846241  0.390085  0.993611
 0.898843  0.774898  0.166567  -0.38555
  0.86173 -0.286951 -0.514705 0.0611899
-0.896957  -0.78822 -0.644434  0.419005
B1_vec: -0.321707  0.728709  0.860527 0.0241207
-0.931782 -0.581732  0.871519 -0.604223
-0.467583  0.833312 -0.386847  0.821305
-0.665156 -0.522945  0.259179   -0.3304
A:  0.964414  0.846241  0.390085  0.993611
 0.898843  0.774898  0.166567  -0.38555
  0.86173 -0.286951 -0.514705 0.0611899
-0.896957  -0.78822 -0.644434  0.419005
B: -0.321707  0.728709  0.860527 0.0241207
-0.931782 -0.581732  0.871519 -0.604223
-0.467583  0.833312 -0.386847  0.821305
-0.665156 -0.522945  0.259179   -0.3304
MeanA:   1.25306  0.478006 -0.166341  0.703709
 0.767439  0.942536  0.419877 -0.253574
 0.219445  0.532435   0.72344  0.706271
 -1.08687 -0.545941 -0.278335  0.609744
MeanB:  -0.283397   0.767624   0.797547  0.0875792
-0.0243491   0.640805  -0.734126   0.793441
  -1.30756  -0.298675    1.09957   -0.47235
 -0.144949   0.113973  -0.637086   0.493252
SigA:  1.53251e-34 -1.35675e-32  9.08441e-33 -1.13558e-33 -1.56173e-34 -2.20257e-33
-1.35675e-32  1.20115e-30 -8.04254e-31  1.00535e-31  1.38262e-32  1.94996e-31
 9.08441e-33 -8.04254e-31  5.38505e-31  -6.7315e-32 -9.25759e-33 -1.30563e-31
-1.13558e-33  1.00535e-31  -6.7315e-32  8.41461e-33  1.15723e-33  1.63209e-32
-1.56173e-34  1.38262e-32 -9.25759e-33  1.15723e-33   1.5915e-34  2.24455e-33
-2.20257e-33  1.94996e-31 -1.30563e-31  1.63209e-32  2.24455e-33  3.16558e-32
SigB:  6.06757e-33  5.67447e-33  -4.6562e-33  5.05973e-32 -2.48815e-32  2.40123e-32
 5.67447e-33  5.30683e-33 -4.35454e-33  4.73192e-32 -2.32695e-32  2.24566e-32
 -4.6562e-33 -4.35454e-33  3.57313e-33 -3.88279e-32  1.90938e-32 -1.84268e-32
 5.05973e-32  4.73192e-32 -3.88279e-32  4.21929e-31 -2.07486e-31  2.00238e-31
-2.48815e-32 -2.32695e-32  1.90938e-32 -2.07486e-31  1.02032e-31  -9.8468e-32
 2.40123e-32  2.24566e-32 -1.84268e-32  2.00238e-31  -9.8468e-32  9.50282e-32
 SigA: 3.10504e+231 -1.35675e-32  9.08441e-33
3.10504e+231         -0.1 -8.04254e-31
 9.08441e-33 -8.04254e-31         -0.1
SigB:  5.05973e-32 -2.48815e-32  2.40123e-32
 4.73192e-32 -2.32695e-32  2.24566e-32
-3.88279e-32  1.90938e-32 -1.84268e-32
Z_g:            1            0            0 -4.40361e-31
           0    -0.995286    0.0969848  1.01371e-30
           0   -0.0969848    -0.995286  1.51396e-30
           0            0            0            0
Z_g.size(): 8
Z_g[0].size(): 16
s_Z: 0
Z_index: 0
Z_g:
           1            0            0 -4.40361e-31
           0    -0.995286    0.0969848  1.01371e-30
           0   -0.0969848    -0.995286  1.51396e-30
           0            0            0            0
0
          -1            0           -0 -1.20339e-30
           0     0.289436    -0.957197  1.01371e-30
           0    -0.957197    -0.289436  1.51396e-30
           0            0            0            0
1
          -1            0           -0 -1.20339e-30
           0    -0.289436     0.957197  3.32594e-31
           0     0.957197     0.289436  4.96727e-31
           0            0            0            0
2
           1            0            0 -4.40361e-31
           0     0.995286   -0.0969848  3.32594e-31
           0    0.0969848     0.995286  4.96727e-31
           0            0            0            0
3
          -1            0            0  4.40361e-31
           0     0.995286   -0.0969848 -1.01371e-30
          -0    0.0969848     0.995286 -1.51396e-30
           0            0            0            0
4
           1            0            0  1.20339e-30
           0    -0.289436     0.957197 -1.01371e-30
           0     0.957197     0.289436 -1.51396e-30
           0            0            0            0
5
           1            0            0  1.20339e-30
           0     0.289436    -0.957197 -3.32594e-31
           0    -0.957197    -0.289436 -4.96727e-31
           0            0            0            0
6
          -1            0           -0  4.40361e-31
           0    -0.995286    0.0969848 -3.32594e-31
           0   -0.0969848    -0.995286 -4.96727e-31
           0            0            0            0
7
Z_final:
works till here - solve for Z? - YES
Build successful?


sum_se:
  3.44183   5.53643  -1.45761     7.049
  3.11675  -4.72253   4.30084  -6.86802
   3.4141   6.67006  -2.64389   7.91566
 -6.61389   -1.8921  -2.74761 0.0761945
diff_se:
 1.15986e-15 -1.16146e-15 -1.41617e-16  3.45586e-16
 6.20089e-15  1.41785e-14  1.90024e-15 -5.58323e-15
-8.61072e-16 -3.37249e-15  4.21548e-16  5.60672e-15
-6.83754e-15 -1.04367e-14  4.56729e-17 -1.54903e-15
X[0]:
 0.964414  0.846241  0.390085  0.993611
 0.898843  0.774898  0.166567  -0.38555
  0.86173 -0.286951 -0.514705 0.0611899
-0.896957  -0.78822 -0.644434  0.419005
Mean.inverse():
 0.388696 -0.965176  0.223105  -1.10841
 0.118085   2.36089 -0.704028   1.66102
 -0.68097  -1.26509   1.42491  -1.39068
 0.487732 -0.184069  0.417766  0.516687
Mean.inverse() * X[0]:
 0.693772  0.390666  0.590321  0.307561
 0.139409  0.822151 -0.268741 -0.140016
  0.68141 -0.869298 -0.313565 -0.684375
 0.201482 -0.257037   -0.3884  0.797641
(Mean.inverse() * X[0]).log():
 1.23795e-17 -5.38176e-17  9.35911e-17 -9.17312e-17
-1.09597e-15 -7.67207e-16 -3.54554e-16 -1.26155e-17
 7.33829e-16  2.33342e-17  2.66597e-16 -1.77921e-16
  1.0649e-16  7.35119e-16 -1.81247e-16  1.55163e-16
diff_se:
 1.23795e-17 -5.38176e-17  9.35911e-17 -9.17312e-17
-1.09597e-15 -7.67207e-16 -3.54554e-16 -1.26155e-17
 7.33829e-16  2.33342e-17  2.66597e-16 -1.77921e-16
  1.0649e-16  7.35119e-16 -1.81247e-16  1.55163e-16
*/
