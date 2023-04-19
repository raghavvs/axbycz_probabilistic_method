/*
DESCRIPTION:

The given code defines a function named Wang_AXBYCZ, which computes the transformation matrices
X, Y, and Z that relate three sets of rigid transformations A, B, and C using Wang's AXB=YCZ
algorithm. The function returns the computed matrices X, Y, Z, and the number of iterations
performed.

The input parameters are:
A, B, and C: vectors of Eigen::Matrix4d objects, each representing a rigid transformation
(rotation and translation).
Xact, Yact, and Zact: Eigen::Matrix4d objects, representing the initial approximations of
the transformation matrices X, Y, and Z.
The function starts by extracting the rotation and translation components from the input
matrices A, B, and C into separate vectors of Eigen::Matrix3d and Eigen::Vector3d objects.
The initial estimates for the rotation matrices RX_init, RY_init, and RZ_init are computed
using skewExp and so3Vec functions.
A while loop is used to perform an iterative process that refines the estimates of the
rotation matrices RX_init, RY_init, and RZ_init. The loop stops when the change in the
estimates (delR) is below a given threshold or the maximum number of iterations (500) is
reached.
Inside the loop, the function constructs the Jacobian matrix (F) and the residual vector
(q) for the current estimates of the rotation matrices. Then, it computes the update (delR)
by solving a least-squares problem using the LDLT decomposition method.
The rotation matrices RX_init, RY_init, and RZ_init are updated by applying the exponential
map (skewExp) using the computed delR values.
After the loop, the function constructs another Jacobian matrix (J) and residual vector (p)
for the translation components. It then solves another least-squares problem using the LDLT
decomposition method to compute the translation vectors tX, tY, and tZ.
Finally, the function constructs the final transformation matrices X, Y, and Z by combining
the computed rotation and translation components. It returns a tuple containing the matrices
X, Y, Z, and the number of iterations performed (n_step).

Input:
A1, B1, C1, A2, B2, C2: Matrices - dim 4x4
opt: bool
nstd1, nst2: standard deviation
Output:
X_final, Y_final, Z_final: Matrices - dim 4x4

In the case of two robotic arms:
 A - robot 1's base to end effector transformation (forward kinematics)
 B - camera to calibration target transformation
 C - robot 2's base to end effector transformation (forward kinematics)
 X - end effector of robot 1 to camera transformation
 Y - robot 1's base to robot 2's base transformation
 Z - end effector of robot 2 to calibration target transformation
*/

#include <eigen3/Eigen/Dense>
#include <vector>
#include <iostream>
#include "so3Vec.h"
#include "skewLog.h"
#include "skewExp.h"

std::tuple<Eigen::Matrix4d, Eigen::Matrix4d, Eigen::Matrix4d, int> Wang_AXBYCZ(const std::vector<Eigen::Matrix4d>& A, const std::vector<Eigen::Matrix4d>& B, const std::vector<Eigen::Matrix4d>& C, const Eigen::Matrix4d& Xact, const Eigen::Matrix4d& Yact, const Eigen::Matrix4d& Zact) {
    int num = A.size();
    std::vector<Eigen::Matrix3d> RA(num), RB(num), RC(num);
    std::vector<Eigen::Vector3d> TA(num), TB(num), TC(num);

    for (int i = 0; i < num; ++i) {
        RA[i] = A[i].block<3, 3>(0, 0);
        RB[i] = B[i].block<3, 3>(0, 0);
        RC[i] = C[i].block<3, 3>(0, 0);
        TA[i] = A[i].block<3, 1>(0, 3);
        TB[i] = B[i].block<3, 1>(0, 3);
        TC[i] = C[i].block<3, 1>(0, 3);
    }

    double e = M_PI / 5;
    Eigen::Matrix3d RX_init = skewExp(so3Vec(Eigen::Vector3d::Constant(e))) * Xact.block<3, 3>(0, 0);
    Eigen::Matrix3d RZ_init = skewExp(so3Vec(Eigen::Vector3d::Constant(e))) * Zact.block<3, 3>(0, 0);
    Eigen::Matrix3d RY_init = RA[0] * RX_init * RB[0] / RZ_init / RC[0];

    Eigen::VectorXd delR = 10000 * Eigen::VectorXd::Ones(9);
    int n_step = 0;
    const double threshold = 0.01;

    while (delR.norm() > threshold && n_step < 500) {
        Eigen::VectorXd q = Eigen::VectorXd::Zero(num * 9);
        Eigen::MatrixXd F = Eigen::MatrixXd::Zero(num * 9, 9);

        for (int i = 0; i < num; ++i) {
            Eigen::Matrix3d tmp1 = RX_init * RB[i];
            Eigen::Matrix3d tmp2 = RY_init * RC[i] * RZ_init;
            Eigen::Matrix3d qq = -RA[i] * tmp1 + tmp2;
            q.segment<9>(i * 9) = Eigen::Map<Eigen::VectorXd>(qq.data(), 9);

            Eigen::Matrix3d F11 = -RA[i] * so3Vec(tmp1.col(0));
            Eigen::Matrix3d F21 = -RA[i] * so3Vec(tmp1.col(1));
            Eigen::Matrix3d F31 = -RA[i] * so3Vec(tmp1.col(2));
            Eigen::Matrix3d F12 = so3Vec(tmp2.col(0));
            Eigen::Matrix3d F22 = so3Vec(tmp2.col(1));
            Eigen::Matrix3d F32 = so3Vec(tmp2.col(2));
            Eigen::Matrix3d F13 = RY_init * RC[i] * so3Vec(RZ_init.col(0));
            Eigen::Matrix3d F23 = RY_init * RC[i] * so3Vec(RZ_init.col(1));
            Eigen::Matrix3d F33 = RY_init * RC[i] * so3Vec(RZ_init.col(2));

            F.block<3, 3>(i * 9, 0) = F11;
            F.block<3, 3>(i * 9, 3) = F12;
            F.block<3, 3>(i * 9, 6) = F13;
            F.block<3, 3>(i * 9 + 3, 0) = F21;
            F.block<3, 3>(i * 9 + 3, 3) = F22;
            F.block<3, 3>(i * 9 + 3, 6) = F23;
            F.block<3, 3>(i * 9 + 6, 0) = F31;
            F.block<3, 3>(i * 9 + 6, 3) = F32;
            F.block<3, 3>(i * 9 + 6, 6) = F33;
        }

        delR = (F.transpose() * F).ldlt().solve(F.transpose() * q);

        double thetaX = delR.segment<3>(0).norm();
        RX_init = skewExp(delR.segment<3>(0) / thetaX, thetaX) * RX_init;

        double thetaY = delR.segment<3>(3).norm();
        RY_init = skewExp(delR.segment<3>(3) / thetaY, thetaY) * RY_init;

        double thetaZ = delR.segment<3>(6).norm();
        RZ_init = skewExp(delR.segment<3>(6) / thetaZ, thetaZ) * RZ_init;

        n_step++;
    }

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3 * num, 9);
    Eigen::VectorXd p = Eigen::VectorXd::Zero(3 * num);

    for (int i = 0; i < num; ++i) {
        J.block<3, 3>(i * 3, 0) = RA[i];
        J.block<3, 3>(i * 3, 3) = -Eigen::Matrix3d::Identity();
        J.block<3, 3>(i * 3, 6) = -RY_init * RC[i];
        p.segment<3>(i * 3) = -TA[i] - RA[i] * RX_init * TB[i] + RY_init * TC[i];
    }

    Eigen::VectorXd translation = (J.transpose() * J).ldlt().solve(J.transpose() * p);
    Eigen::Vector3d tX = translation.segment<3>(0);
    Eigen::Vector3d tY = translation.segment<3>(3);
    Eigen::Vector3d tZ = translation.segment<3>(6);

    Eigen::Matrix4d X, Y, Z;
    X.setIdentity();
    Y.setIdentity();
    Z.setIdentity();
    X.block<3, 3>(0, 0) = RX_init;
    Y.block<3, 3>(0, 0) = RY_init;
    Z.block<3, 3>(0, 0) = RZ_init;
    X.block<3, 1>(0, 3) = tX;
    Y.block<3, 1>(0, 3) = tY;
    Z.block<3, 1>(0, 3) = tZ;

    return std::make_tuple(X, Y, Z, n_step);
}

int main() {
    // Define A, B, C, Xact, Yact, and Zact matrices, and populate them with appropriate data
    std::vector<Eigen::Matrix4d> A, B, C;
    Eigen::Matrix4d Xact, Yact, Zact;

    // Populate A, B, C, Xact, Yact, and Zact with data

    auto [X, Y, Z, n_step] = Wang_AXBYCZ(A, B, C, Xact, Yact, Zact);

    std::cout << "X: \n" << X << std::endl;
    std::cout << "Y: \n" << Y << std::endl;
    std::cout << "Z: \n" << Z << std::endl;
    std::cout << "Number of steps: " << n_step << std::endl;

    return 0;
}