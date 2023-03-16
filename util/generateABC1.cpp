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

std::tuple<std::vector<Eigen::Matrix4d>, std::vector<Eigen::Matrix4d>, std::vector<Eigen::Matrix4d>>
    generateABC(int length, int optFix, int optPDF, const Eigen::VectorXd& M, const Eigen::MatrixXd& Sig,
                const Eigen::Matrix4d& X, const Eigen::Matrix4d& Y, const Eigen::Matrix4d& Z)
{
    std::vector<Eigen::Matrix4d> A(length), B(length), C(length);

    Eigen::VectorXd a = Eigen::VectorXd::Random(6).normalized();
    Eigen::Matrix4d A_initial = (Eigen::Matrix4d(se3Vec(a))).exp();

    Eigen::VectorXd b = Eigen::VectorXd::Random(6).normalized();
    Eigen::Matrix4d B_initial = (Eigen::Matrix4d(se3Vec(b))).exp();

    Eigen::VectorXd c = Eigen::VectorXd::Random(6).normalized();
    Eigen::Matrix4d C_initial = (Eigen::Matrix4d(se3Vec(c))).exp();

    if (optFix == 1) { // Fix A, randomize B and C - This can be applied to both serial-parallel and dual-robot arm calibrations
        for (int m = 0; m < length; m++) {
            if (optPDF == 1) {
                Eigen::VectorXd randVec = mvg(M, Sig, 1).first;
                // Update B matrix with random noise
                B[m] = (Eigen::MatrixXd(se3Vec(randVec)).exp() * B_initial);
            }
            // Compute C matrix from A,B,X,Y,Z matrices
            C[m] = Y.inverse() * (A_initial * X * B[m] * Z.inverse());
            // Assign fixed value to A matrix
            A[m] = A_initial;
        }
    }

    // Return a tuple of vectors containing the output matrices
    return std::make_tuple(A,B,C);
}

int main() {

    int length = 2;
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

    // Call the function and get the output matrices in a tuple
    auto [A,B,C] = generateABC(length, optFix, optPDF, M, Sig, X, Y, Z);

    for(int i = 0; i < length; i++) {
        std::cout << "A: \n" << A[i] << std::endl;
        std::cout << "B: \n" << B[i] << std::endl;
        std::cout << "C: \n" << C[i] << std::endl;
    }

    return 0;
}