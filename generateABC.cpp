#include <iostream>
#include <cmath
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <random>

using namespace Eigen;

void generateABC(int length, int optFix, int optPDF, VectorXd M, MatrixXd Sig, Matrix4d X, Matrix4d Y, Matrix4d Z, Matrix4d& A, Matrix4d& B, Matrix4d& C) {

    int len = length;
    int dataGenMode = 3;
    double qz1[6] = {M_PI/6, M_PI/3, M_PI/4, M_PI/4, -M_PI/4, 0};
    double qz2[6] = {M_PI/3, M_PI/4, M_PI/3, -M_PI/4, M_PI/4, 0};
    double qz3[6] = {M_PI/4, M_PI/3, M_PI/3, M_PI/6, -M_PI/4, 0};
    Matrix<double, 6, 1> a, b, c;
    Matrix4d A_initial, B_initial, C_initial;

    if (dataGenMode == 1) {

        // Instantiate puma object (time-consuming)
        std::cout << "Error: Puma object not instantiated. Exiting program." << std::endl;
        exit(1);

    } else if (dataGenMode == 2) {

        A << 0.2294, -0.1951, -0.9536, -0.1038,
             0.7098,  0.7039,  0.0268, -0.2332,
             0.6660, -0.6830,  0.3000,  0.2818,
             0.0,     0.0,     0.0,     1.0;
        B << 0.0268, -0.7039, -0.7098,  0.0714,
            -0.9536,  0.1951, -0.2294, -0.1764,
             0.3000,  0.6830, -0.6660,  0.2132,
             0.0,     0.0,     0.0,     1.0;
        C << -0.0335, -0.4356, -0.8995, -0.0128,
             0.4665,  0.7891, -0.3995, -0.2250,
             0.8839, -0.4330,  0.1768,  0.1756,
             0.0,     0.0,     0.0,     1.0;

    } else if (dataGenMode == 3) {

        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0, 1.0);
        auto randn = [&] { return distribution(generator); };

        a << randn(), randn(), randn(), randn(), randn(), randn();
        a.normalize();
        A_initial = Matrix4d::Identity();
        A_initial.topLeftCorner<3, 3>() = AngleAxisd(a.norm(), a.normalized()).toRotationMatrix();

        b << randn(), randn(), randn(), randn(), randn(), randn();
        b.normalize();
        B_initial = Matrix4d::Identity();
        B_initial.topLeftCorner<3, 3>() = AngleAxisd(b.norm(), b.normalized()).toRotationMatrix();

        c << randn(), randn(), randn(), randn(), randn(), randn();
        c.normalize();
        C_initial = Matrix4d::Identity();
        C_initial.topLeftCorner<3, 3>() = AngleAxisd(c.norm(), c.normalized()).toRotationMatrix();
    }

       
