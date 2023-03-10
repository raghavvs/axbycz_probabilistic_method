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
*/

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include <random>
#include <se3Vec.h>
#include <mvg.h>
#include <sensorNoise.h>
#include <fKine.h>


void generateABC(int length, int optFix, int optPDF, Eigen::VectorXd M, Eigen::MatrixXd Sig, 
                Eigen::Matrix4d X, Eigen::Matrix4d Y, Eigen::Matrix4d Z, Eigen::Matrix4d& A, 
                Eigen::Matrix4d& B, Eigen::Matrix4d& C) {
    int len = length;
    int dataGenMode = 3;
    Eigen::Matrix4d A_initial, B_initial, C_initial;
    const Eigen::VectorXd qz1 = {M_PI/6, M_PI/3, M_PI/4, M_PI/4, -M_PI/4, 0};
    const Eigen::VectorXd qz2 = {M_PI/3, M_PI/4, M_PI/3, -M_PI/4, M_PI/4, 0};
    const Eigen::VectorXd qz3 = {M_PI/4, M_PI/3, M_PI/3, M_PI/6, -M_PI/4, 0};
    Eigen::Matrix<double, 6, 1> a, b, c;

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

        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0, 1.0);
        auto randn = [&] { return distribution(generator); };

        a << randn(), randn(), randn(), randn(), randn(), randn();
        a.normalize();
        A_initial = Eigen::Matrix4d::Identity();
        A_initial.topLeftCorner<3, 3>() = Eigen::AngleAxisd(a.norm(), a.normalized()).toRotationMatrix();

        b << randn(), randn(), randn(), randn(), randn(), randn();
        b.normalize();
        B_initial = Matrix4d::Identity();
        B_initial.topLeftCorner<3, 3>() = Eigen::AngleAxisd(b.norm(), b.normalized()).toRotationMatrix();

        c << randn(), randn(), randn(), randn(), randn(), randn();
        c.normalize();
        C_initial = Matrix4d::Identity();
        C_initial.topLeftCorner<3, 3>() = Eigen::AngleAxisd(c.norm(), c.normalized()).toRotationMatrix();
    }
    if (optFix == 1) { // Fix A, randomize B and C
    // This can be applied to both serial-parallel and dual-robot arm
    // calibrations

    Eigen::MatrixXd B_initial(4, 4);
    B_initial << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1;

    Eigen::MatrixXd C_initial(4, 4);
    C_initial << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1;

    Eigen::MatrixXd X(4, 4);
    Eigen::MatrixXd Y(4, 4);
    Eigen::MatrixXd Z(4, 4);

    Eigen::MatrixXd B(4, 4, len);
    Eigen::MatrixXd C(4, 4, len);
    Eigen::MatrixXd A(4, 4, len);

    for (int m = 0; m < len; m++) {

        if (optPDF == 1) {
            Eigen::VectorXd randn_vec = mvg(M, Sig, 1);
            Eigen::MatrixXd expm_vec = se3Vec(randn_vec);
            B.block<4, 4>(0, 0, 4, 4) = (expm_vec).exp() * B_initial;
        }
        else if (optPDF == 2) {
            Eigen::VectorXd randn_vec = mvg(M, Sig, 1);
            Eigen::MatrixXd expm_vec = se3Vec(randn_vec);
            B.block<4, 4>(0, 0, 4, 4) = B_initial * (expm_vec).exp();
        }
        else if (optPDF == 3) {
            Eigen::VectorXd gmean(6);
            gmean << 0, 0, 0, 0, 0, 0;
            double sigma = Sig(0, 0);
            B.block<4, 4>(0, 0, 4, 4) = sensorNoise(B_initial, gmean, sigma, 1);
        }

        C.block<4, 4>(0, 0, 4, 4) = Y.inverse() * (A_initial * X * B.block<4, 4>(0, 0, 4, 4) / Z);
        A.block<4, 4>(0, 0, 4, 4) = A_initial;
        }
    }

    if (optFix == 2) { // Fix B, randomize A and C
    // This can be applied to both serial-parallel and dual-robot arm calibrations
    Eigen::MatrixXd A(4, 4, len);
    Eigen::MatrixXd C(4, 4, len);

    for (int m = 0; m < len; m++) {
        if (optPDF == 1) {
            Eigen::VectorXd random_vector = mvg(M, Sig, 1);
            Eigen::MatrixXd matrix = se3Vec(random_vector);
            A.slice(m) = matrix.exp() * A_initial;
        } else if (optPDF == 2) {
            Eigen::VectorXd random_vector = mvg(M, Sig, 1);
            Eigen::MatrixXd matrix = se3Vec(random_vector);
            A.slice(m) = A_initial * matrix.exp();
        } else if (optPDF == 3) {
            Eigen::VectorXd gmean = Eigen::VectorXd::Zero(6);
            A.slice(m) = sensorNoise(A_initial, gmean, Sig(0), 1);
        }
        
        C.slice(m) = Y.llt().solve(A.slice(m) * X * B_initial) / Z;
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
            B(m) = sensorNoise(B_initial, gmean, Sig(0), 1);
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