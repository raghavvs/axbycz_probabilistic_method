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

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> dist(0.0, 1.0);

        VectorXd a(6);
        a << dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen);
        a = a / a.norm();
        A = Matrix4d::Identity();
        A.block<3,3>(0,0) = AngleAxisd(a.tail(3).norm(), a.tail(3).normalized()).toRotationMatrix();
        A.block<3,1>(0,3) = a.head(3);

        VectorXd b(6);
        b << dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen);
       
