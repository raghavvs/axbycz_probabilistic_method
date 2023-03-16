/*
DESCRIPTION:

The program generates three homogeneous transformation matrices 
A, B, and C, each of which is represented as a 4x4 matrix. The 
matrices are generated in a variety of ways depending on the values 
of several input parameters. The matrices can be used to represent 
the positions and orientations of robot arms, and the program includes 
functions for fixing and randomizing the position and orientation 
of the matrices. The program relies on several external libraries, 
including Eigen, which is a library for linear algebra in C++, and 
mvg and sensorNoise, which are libraries for computer vision and 
sensor noise modeling, respectively.

Data generation for AXB = YCZ problem
Input:
       length: number of generated data pairs
       optFix: option for fixing different data streams
       optPDF: option for generating data using different distributions
       M:      mean of perturbance in lie algebra
       Sig:    covariance of perturbance in lie algebra
       X, Y, Z: ground truths
 Output:
       A, B, C: 4 x 4 x length or 4 x 4
                noise-free data streams with correspondence
*/

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include <random>
#include "mvg.h"
#include "sensorNoise.h"
//#include "se3Vec.h"
#include "fKine.h"

void generateABC(int length, int optFix, int optPDF, const Eigen::VectorXd M, const Eigen::MatrixXd Sig,
                 const Eigen::Matrix4d X, const Eigen::Matrix4d Y, Eigen::Matrix4d Z, const Eigen::Matrix4d& A,
                 const Eigen::Matrix4d& B, const Eigen::Matrix4d& C)
{
    int len = length;
    Eigen::Matrix4d A_initial, B_initial, C_initial;

    Eigen::VectorXd a = Eigen::VectorXd::Random(6);
    a /= a.norm();
    A_initial = (Eigen::Matrix4d(se3Vec(a))).exp();

    Eigen::VectorXd b = Eigen::VectorXd::Random(6);
    b /= b.norm();
    B_initial = (Eigen::Matrix4d(se3Vec(b))).exp();

    Eigen::VectorXd c = Eigen::VectorXd::Random(6);
    c /= c.norm();
    C_initial = (Eigen::Matrix4d(se3Vec(c))).exp();

    //PART II - Fix a matrix A, B, C - Only using Gaussian noise - optPDF = 1

    if (optFix ==1) { // Fix A, randomize B and C - This can be applied to both serial-parallel and dual-robot arm calibrations
        Eigen::Matrix4d A[len], B[len], C[len];
        for (int m = 0; m < len; m++) {
            if (optPDF == 1) {
                Eigen::Matrix<double, 6, 1> randVec = mvg(M, Sig, 1).first;
                B[m] = (Eigen::Matrix4d(se3Vec(randVec)).exp() * B_initial);
            }

            C[m] = Y.inverse() * (A_initial * X * B[m] * Z.inverse());
            A[m] = A_initial;
        }
    }
}

int main() {
    int length = 10;
    int optFix = 1;
    int optPDF = 1;
    Eigen::VectorXd M(6);
    M << 0, 0, 0, 0, 0, 0;
    Eigen::MatrixXd Sig(6, 6);
    Sig << 0.1, 0, 0, 0, 0, 0,
            0, 0.1, 0, 0, 0, 0,
            0, 0, 0.1, 0, 0, 0,
            0, 0, 0, 0.01, 0, 0,
            0, 0, 0, 0, 0.01, 0,
            0, 0, 0, 0, 0, 0.01;
    Eigen::Matrix4d X = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d Y = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d Z = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d A, B, C = Eigen::Matrix4d::Zero();
    generateABC(length, optFix, optPDF, M, Sig, X, Y, Z, A, B, C);
    std::cout << "A: \n" << A << std::endl;
    std::cout << "B: \n" << B << std::endl;
    std::cout << "C: \n" << C << std::endl;
    return 0;
}
