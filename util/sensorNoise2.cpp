#include <iostream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include <random>
#include "se3Vec.h"
#include "so3Vec.h"

Eigen::MatrixXd sensorNoise(const Eigen::Matrix4d g[], Eigen::MatrixXd gmean, double sd) {
    
    Eigen::Matrix4d g_noise[10];

    std::cout << "g: " << g->rows() << " x " << g->cols() << std::endl;

    Eigen::Vector3d temp = Eigen::Vector3d::Random();
    Eigen::VectorXd noise_old1(6);
    Eigen::VectorXd noise_old2(6);
    Eigen::Matrix4d exp1;
    Eigen::Matrix4d exp2;
    Eigen::Matrix4d prod;

    // Independently from Normal Distribution
    noise_old1.segment(0, 3) = Eigen::Vector3d::Zero();
    noise_old1.segment(3, 3) = sd * Eigen::Vector3d::Random();
    noise_old1 += gmean;

    noise_old2.segment(0, 3) = sd * temp.normalized() * temp.norm();
    noise_old2.segment(3, 3) = Eigen::Vector3d::Zero();
    noise_old2 += gmean;

    exp1 = (se3Vec(noise_old1)).exp();
    exp2 = (se3Vec(noise_old2)).exp();
    prod = exp1 * exp2;

    for (int i = 0; i < 10; i++) {
        //g_noise[i] = g[i] * exp1[i] * exp2[i];
        g_noise[i] = g[i] * prod;
    }

    std::cout << "g_noise: " << g_noise->rows() << " x " << g_noise->cols() << std::endl;

    for (int i = 0; i < 10; i++) {
        std::cout << g_noise[i] << std::endl;
    };

    std::cout << "works" << std::endl;

    Eigen::Matrix4d a = Eigen::Matrix4d::Random();

    return a;
    //return g_noise;
}

int main() {

    // Create an array of 10 matrices of size 4x4
    Eigen::Matrix4d m[10];

    // Generate random matrices and store them in the array
    for (int i = 0; i < 10; i++) {
        m[i] = Eigen::Matrix4d::Random();
    }

    double sd = 0.1;

    // set noise variables
    Eigen::MatrixXd g_mean(6, 1);
    g_mean << 0, 0, 0, 0, 0, 0;

    // apply noise to matrix and print dimensions
    Eigen::MatrixXd m_noise = sensorNoise(m, g_mean, sd);

    return 0;
}