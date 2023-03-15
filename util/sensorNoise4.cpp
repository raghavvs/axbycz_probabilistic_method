#include <iostream>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include "se3Vec.h"
#include "so3Vec.h"

Eigen::Matrix4d* sensorNoise(const Eigen::Matrix4d* g, int len, const Eigen::MatrixXd& gmean, double sd) {

    // Declare g_noise as an array of matrices and allocate memory for it
    Eigen::Matrix4d* g_noise = new Eigen::Matrix4d[len];

    //std::cout << "g: " << g->rows() << " x " << g->cols() << std::endl;

    Eigen::Vector3d temp = Eigen::Vector3d::Random();
    Eigen::VectorXd noise_old1(6);
    Eigen::VectorXd noise_old2(6);

    // Independently from Normal Distribution
    noise_old1.segment(0, 3) = Eigen::Vector3d::Zero();
    noise_old1.segment(3, 3) = sd * Eigen::Vector3d::Random();
    noise_old1 += gmean;

    noise_old2.segment(0, 3) = sd * temp.normalized() * temp.norm();
    noise_old2.segment(3, 3) = Eigen::Vector3d::Zero();
    noise_old2 += gmean;

    Eigen::Matrix4d exp1 = (se3Vec(noise_old1)).exp();
    Eigen::Matrix4d exp2 = (se3Vec(noise_old2)).exp();
    Eigen::Matrix4d prod = exp1 * exp2;

    for (int i = 0; i < len; i++) {
        g_noise[i] = g[i] * prod;
    }

    for (int i = 0; i < len; i++) {
        std::cout << "Matrix: " << i << std::endl;
        std::cout << g_noise[i] << std::endl;
    };

    std::cout << "function works" << std::endl;

    return g_noise;
}

int main() {

    Eigen::Matrix4d m[5];
    int size = sizeof(m) / sizeof(m[0]);

    // Generate random matrices and store them in the array
    for (int i = 0; i < size; i++) {
        m[i] = Eigen::Matrix4d::Random();
    }

    double sd = 0.1;

    // set noise variables
    Eigen::MatrixXd g_mean(6, 1);
    g_mean << 0, 0, 0, 0, 0, 0;

    // apply noise to matrix and print dimensions
    Eigen::Matrix4d* m_noise_new = sensorNoise(m, size, g_mean, sd);

    for (int i = 0; i < size; i++) {
        std::cout << "Matrix: " << i << std::endl;
        std::cout << m_noise_new[i] << std::endl;
    };

    delete[] m_noise_new;
    std::cout << size << std::endl;
    std::cout << "main works" << std::endl;

    return 0;

}