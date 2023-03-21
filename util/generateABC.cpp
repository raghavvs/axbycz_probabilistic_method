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
#include "se3Vec.h"
#include "mvg.h"
#include "sensorNoise.h"
#include "fKine.h"

void generateABC(int length, int optFix, int optPDF, Eigen::VectorXd M, Eigen::MatrixXd Sig, 
                Eigen::Matrix4d X, Eigen::Matrix4d Y, Eigen::Matrix4d Z, Eigen::Matrix4d& A, 
                Eigen::Matrix4d& B, Eigen::Matrix4d& C) 
{
    int len = length;
    int dataGenMode = 3;
    
    Eigen::Matrix4d A_initial, B_initial, C_initial;
    Eigen::VectorXd qz1(6);
    qz1 << M_PI/6, M_PI/3, M_PI/4, M_PI/4, -M_PI/4, 0;
    Eigen::VectorXd qz2(6);
    qz2 << M_PI/3, M_PI/4, M_PI/3, -M_PI/4, M_PI/4, 0;
    Eigen::VectorXd qz3(6);
    qz3 << M_PI/4, M_PI/3, M_PI/3, M_PI/6, -M_PI/4, 0;

    Eigen::Matrix<double, 6, 1> a, b, c;

    //PART I -Generate random matrices - A, B, C

    if (dataGenMode == 1) {
        
        A_initial = fKine(qz1);
        B_initial = fKine(qz2);
        C_initial = fKine(qz3);

    } else if (dataGenMode == 2) {

        A_initial << 0.2294, -0.1951, -0.9536, -0.1038,
            0.7098,  0.7039,  0.0268, -0.2332,
            0.6660, -0.6830,  0.3000,  0.2818,
            0.0,     0.0,     0.0,     1.0;
        B_initial << 0.0268, -0.7039, -0.7098,  0.0714,
            -0.9536,  0.1951, -0.2294, -0.1764,
            0.3000,  0.6830, -0.6660,  0.2132,
            0.0,     0.0,     0.0,     1.0;
        C_initial << -0.0335, -0.4356, -0.8995, -0.0128,
            0.4665,  0.7891, -0.3995, -0.2250,
            0.8839, -0.4330,  0.1768,  0.1756,
            0.0,     0.0,     0.0,     1.0;

    } else if (dataGenMode == 3) {

        Eigen::VectorXd a = Eigen::VectorXd::Random(6);
        a /= a.norm(); // normalize a
        Eigen::Matrix4d A_initial = (Eigen::Matrix4d(se3Vec(a))).exp();

        Eigen::VectorXd b = Eigen::VectorXd::Random(6);
        b /= b.norm(); // normalize b
        Eigen::Matrix4d B_initial = (Eigen::Matrix4d(se3Vec(b))).exp();

        Eigen::VectorXd c = Eigen::VectorXd::Random(6);
        c /= c.norm(); // normalize c
        Eigen::Matrix4d C_initial = (Eigen::Matrix4d(se3Vec(c))).exp();

    //PART II - Fix a matrix A, B, C - Only using Gaussian noise - optPDF = 1

    if (optFix == 1){ // Fix A, randomize B and C - This can be applied to both serial-parallel and dual-robot arm calibrations
        Eigen::Matrix4d A[len], B[len], C[len];
        for (int m = 0; m < len; m++){
            if (optPDF == 1){
                Eigen::Matrix<double, 6, 1> randVec = mvg(M, Sig, 1).first;
                B[m] = (Eigen::Matrix4d(se3Vec(randVec)).exp() * B_initial);
            }
            /*else if (optPDF == 2){
                Eigen::Matrix<double, 6, 1> randVec = mvg(M, Sig, 1);
                B[m] = Eigen::Matrix4d(B_initial * (se3Vec(randVec)).exp());
            }
            else if (optPDF == 3){
                Eigen::Matrix<double, 6, 1> g_mean;
                g_mean << 0, 0, 0, 0, 0, 0;
                double sd = Sig(0); // Assuming Sig is a vector with one value
                std::pair<Eigen::Matrix4d*, double> sensor_output = sensorNoise(B_initial, g_mean, sd, 1);
                B[m] = *(sensor_output.first);
            }*/

            C[m] = Y.inverse() * (A_initial * X * B[m] * Z.inverse());
            A[m] = A_initial;
        }
    } else if(optFix == 2) // Fix B, randomize A and C - This can be applied to both serial-parallel and dual-robot arm calibrations
        Eigen::Matrix4d A[len], B[len], C[len];
        for(int m = 0; m < len; m++) {
            if(optPDF == 1) {
                Eigen::Matrix<double, 6, 1> randVec = mvg(M, Sig, 1).first;
                A[m] = (Eigen::Matrix4d(se3Vec(randVec)).exp() * A_initial);
            /*} else if(optPDF == 2) {
                A[m] = A_initial * (se3Vec(mvg(M, Sig, 1))).exp();
            } else if(optPDF == 3) {
                Eigen::VectorXd gmean(6);
                gmean << 0, 0, 0, 0, 0, 0;
                A[m] = sensorNoise(A_initial, gmean, Sig(0,0), 1);
            }*/

            C[m] = Y.inverse() * (A_initial * X * B[m] * Z.inverse());
            B[m] = B_initial;
        }
    }
    Eigen::MatrixXd B = B_initial;

    if (optFix == 3) { // Fix C, randomize A and B
    // This is only physically achievable on multi-robot hand-eye calibration
    Eigen::MatrixXd A(len, 4, 4);
    Eigen::MatrixXd B(len, 4, 4);
    Eigen::MatrixXd C(len, 4, 4);

    for (int m = 0; m < len; m++) {
        if (optPDF == 1) {
            B(m) = Eigen::Matrix4d::Identity() + se3Vec(mvg(M, Sig, 1));
            B(m) *= B_initial;
        }
        else if (optPDF == 2) {
            B(m) = B_initial * (Eigen::Matrix4d::Identity() + se3Vec(mvg(M, Sig, 1)));
        }
        else if (optPDF == 3) {
            Eigen::VectorXd gmean(6);
            gmean << 0, 0, 0, 0, 0, 0;
            B(m) = sensorNoise(B_initial, len, gmean, Sig(0), 1);
        }

        A(m) = (Y * C_initial * Z / B(m)) / X;
        C(m) = C_initial;
    }

    // Update B with the randomized values
    B_initial = B(len-1);
    }

    else if (optFix == 4) {
        // This is for testing tranditional AXBYCZ solver that demands the
        // correspondence between the data pairs {A_i, B_i, C_i}
        Eigen::MatrixXd A(len, 4, 4);
        Eigen::MatrixXd B(len, 4, 4);
        Eigen::MatrixXd C(len, 4, 4);

        for (int m = 0; m < len; m++) {
            A(m) = Eigen::Matrix4d::Identity() + se3Vec(mvg(M, Sig, 1));
            A(m) *= A_initial;

            C(m) = Eigen::Matrix4d::Identity() + se3Vec(mvg(M, Sig, 1));
            C(m) *= C_initial;

            B(m) = X.lu().solve(A(m).lu().solve(Y * C(m) * Z));
        }
    }
}

/* int main() {
    int length = 10;
    int optFix = 1;
    int optPDF = 3;
    Eigen::VectorXd M(6);
    M << 0, 0, 0, 0, 0, 0;
    Eigen::MatrixXd Sig(6, 6);
    Sig << Eigen::MatrixXd::Identity(6, 6);
    Eigen::Matrix4d X = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d Y = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d Z = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d A, B, C;

    generateABC(length, optFix, optPDF, M, Sig, X, Y, Z, A, B, C);

    for (int i = 0; i < length; ++i) {
        std::cout << "A_" << i << ":\n" << A.col(i) << "\n\n";
        std::cout << "B_" << i << ":\n" << B.col(i) << "\n\n";
        std::cout << "C_" << i << ":\n" << C.col(i) << "\n\n";
    }

    return 0;
} */