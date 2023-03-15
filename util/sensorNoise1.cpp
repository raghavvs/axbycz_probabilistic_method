#include <iostream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include <random>
#include "se3Vec.h"
#include "so3Vec.h"

Eigen::MatrixXd sensorNoise(Eigen::MatrixXd g, Eigen::MatrixXd gmean, double sd) {
    Eigen::MatrixXd g_noise(g.rows(), g.cols());

    // Independently from Normal Distribution
    Eigen::Vector3d temp = Eigen::Vector3d::Random();
    Eigen::VectorXd noise_old1(6);
    Eigen::VectorXd noise_old2(6);

    noise_old1.segment(0, 3) = Eigen::Vector3d::Zero();
    noise_old1.segment(3, 3) = sd * Eigen::Vector3d::Random();
    noise_old1 += gmean;

    noise_old2.segment(0, 3) = sd * temp.normalized() * temp.norm();
    noise_old2.segment(3, 3) = Eigen::Vector3d::Zero();
    noise_old2 += gmean;

    std::cout << "g: " << g.rows() << " x " << g.cols() << std::endl;
    std::cout << "noise_old1: " << noise_old1.size() << std::endl;

    Eigen::MatrixXd exp1 = (Eigen::MatrixXd((se3Vec(noise_old1))).exp());
    std::cout << "exp1: " << exp1.rows() << " x " << exp1.cols() << std::endl;

    Eigen::MatrixXd exp2 = (Eigen::MatrixXd((se3Vec(noise_old2))).exp());
    std::cout << "exp2: " << exp2.rows() << " x " << exp2.cols() << std::endl;

    g_noise = g * exp1 * exp2;
    std::cout << "g_noise: " << g_noise.rows() << " x " << g_noise.cols() << std::endl;

    return g_noise;
}

int main() {
    // create a matrix and a vector for testing
    Eigen::Matrix4d m;
    m << 0, -3, 2, 4,
        3, 0, -1, 5,
        -2, 1, 0, 6,
        0, 0, 0, 0;

    // set noise variables
    Eigen::MatrixXd g_mean(6, 1);
    g_mean << 0, 0, 0, 0, 0, 0;

    double sd = 0.1;

    // apply noise to matrix and print dimensions
    Eigen::MatrixXd m_noise = sensorNoise(m, g_mean, sd);
    std::cout << "m_noise: " << m_noise.rows() << " x " << m_noise.cols() << std::endl;
    std::cout << "m_noise: " << std::endl << m_noise << std::endl;
    std::cout << "size: " <<  m_noise.size() << std::endl;

    return 0;
}