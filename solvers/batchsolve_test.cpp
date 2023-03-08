#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
#include <meanCov.h>

using namespace Eigen;

std::vector<MatrixXd> eig(MatrixXd A) {
    EigenSolver<MatrixXd> es(A);
    std::vector<MatrixXd> result(2);
    result[0] = es.eigenvalues().real();
    result[1] = es.eigenvectors().real();
    return result;
}

void batchSolveXY(MatrixXd A, MatrixXd B, bool opt, double nstd_A, double nstd_B, MatrixXd& X, MatrixXd& Y) {

    Tensor<double, 3> X_tensor(4, 4, 8);
    X_tensor.setZero();
    Tensor<double, 3> Y_tensor(4, 4, 8);
    Y_tensor.setZero();

    Matrix4d X_candidate, Y_candidate;

    int a1 = A.rows(), a2 = A.cols(), a3 = A.size() / (a1 * a2);
    Matrix3d A_mex = Map<Matrix3d>(A.data(), a1, a2 * a3);
    Matrix3d B_mex = Map<Matrix3d>(B.data(), a1, a2 * a3);

    MatrixXd MeanA = MatrixXd::Zero(4, 1);
    MatrixXd SigA = MatrixXd::Zero(6, 6);
    MatrixXd MeanB = MatrixXd::Zero(4, 1);
    MatrixXd SigB = MatrixXd::Zero(6, 6);

    meanCov(A_mex, MeanA, SigA);
    meanCov(B_mex, MeanB, SigB);

    if (opt) {
        SigA -= nstd_A * MatrixXd::Identity(6, 6);
        SigB -= nstd_B * MatrixXd::Identity(6, 6);
    }

    MatrixXd VA, VB;
    VA.setZero(3, 3);
    VB.setZero(3, 3);

    std::vector<MatrixXd> VA_eig, VB_eig;
    VA_eig = eig(SigA.block<3, 3>(0, 0));
    VB_eig = eig(SigB.block<3, 3>(0, 0));
    VA = VA_eig[1];
    VB = VB_eig[1];

    MatrixXd Q1, Q2, Q3, Q4;
    Q1 = MatrixXd::Identity(3, 3);
    Q2 = (MatrixXd(3, 3) << -1, 0, 0, 0, -1, 0, 0, 0, 1).finished();
    Q3 = (MatrixXd(3, 3) << -1, 0, 0, 0, 1, 0, 0, 0, -1).finished();
    Q4 = (MatrixXd(3, 3) << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished();

    MatrixXd Rx_solved(3, 3, 8);

    // There are 8 possibilities of Rx
    Rx_solved.block(0, 0, 3, 3) = VA * Q1 * VB.transpose();
    Rx_solved.block(0, 3, 3, 3) = VA * Q2 * VB.transpose();
    Rx_solved.block(3, 0, 3, 3) = -VB * Q1.transpose() * VA.transpose();
    Rx_solved.block(3, 3, 3, 3) = -VB * Q2.transpose() * VA.transpose();
    Rx_solved.block(6, 0, 2, 3) = Eigen::Matrix<double, 2, 3>::Zero();
    Rx_solved.block(6, 3, 2, 3) = Eigen::Matrix<double, 2, 3>::Zero();
    Rx_solved(8, 8) = 1.0;

    // Compute t using equation (22)
    Eigen::Matrix<double, 3, 1> t_solved;
    t_solved = VAi - Rx_solved.block(0, 0, 3, 3) * VB.transpose() * VB1;

    // Assemble the transformation matrix
    T_solved.block(0, 0, 3, 3) = Rx_solved;
    T_solved.block(0, 3, 3, 1) = t_solved;
    T_solved.block(3, 0, 1, 4) = Eigen::Matrix<double, 1, 4>::Zero();
    T_solved(3, 3) = 1.0;

    return T_solved;
}

int main()
{
    // Define two random 3D point clouds
    Eigen::Matrix<double, 3, Eigen::Dynamic> cloud1(3, 6);
    Eigen::Matrix<double, 3, Eigen::Dynamic> cloud2(3, 6);
    cloud1 << 3.0, 2.0, 4.0, 6.0, 8.0, 10.0,
            1.0, 5.0, 4.0, 3.0, 9.0, 7.0,
            0.0, 2.0, 1.0, 3.0, 6.0, 8.0;

    cloud2 << 3.1, 2.2, 3.8, 5.9, 8.1, 9.8,
            1.1, 5.2, 4.1, 2.9, 9.2, 7.1,
            0.1, 1.9, 1.2, 3.1, 6.1, 7.9;

    // Estimate the rigid body transformation between the two point clouds
    Eigen::Matrix4d T = estimateRigidTransformation(cloud1, cloud2);

    std::cout << "Rigid body transformation matrix:" << std::endl;
    std::cout << T << std::endl;

    return 0;
}
