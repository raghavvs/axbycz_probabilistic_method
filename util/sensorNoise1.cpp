#include <iostream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include <random>
#include <se3Vec.h>
#include <so3Vec.h>

using namespace Eigen;

MatrixXd sensorNoise(MatrixXd g, MatrixXd gmean, double std, int model) {
    MatrixXd g_noise(g.rows(), g.cols());

    switch (model) {
        case 1:
            // Independently from Normal Distribution
            for (int i = 0; i < g.cols(); i++) {
                Vector3d temp = Vector3d::Random();
                VectorXd noise_old1(6);
                VectorXd noise_old2(6);

                noise_old1.segment(0, 3) = Vector3d::Zero();
                noise_old1.segment(3, 3) = std * Vector3d::Random();
                noise_old1 += gmean;

                noise_old2.segment(0, 3) = std * temp.normalized() * temp.norm();
                noise_old2.segment(3, 3) = Vector3d::Zero();
                noise_old2 += gmean;

                std::cout << "g: " << g.rows() << " x " << g.cols() << std::endl;
                std::cout << "noise_old1: " << noise_old1.size() << std::endl;
                
                MatrixXd exp1 = (MatrixXd((se3Vec(noise_old1))).exp());
                std::cout << "exp1: " << exp1.rows() << " x " << exp1.cols() << std::endl;

                MatrixXd exp2 = (MatrixXd((se3Vec(noise_old2))).exp());
                std::cout << "exp2: " << exp2.rows() << " x " << exp2.cols() << std::endl;

                g_noise.col(i) = g.col(i) * exp1 * exp2.topRows(g.rows());
                std::cout << "g_noise: " << g_noise.rows() << " x " << g_noise.cols() << std::endl;
            }
            break;
        default:
            std::cerr << "Invalid noise model specified" << std::endl;
            break;
    }

    return g_noise;
}

int main() {
    // create a matrix and a vector for testing
    Eigen::MatrixXd m(3, 2);
    m << 1, 2,
         3, 4,
         5, 6;

    Eigen::VectorXd v(3);
    v << 7, 8, 9;

    // set noise variables
    Eigen::MatrixXd g_mean(4, 1);
    g_mean << 0, 0, 0, 0;

    double std_dev = 0.1;
    int model = 1;

    // apply noise to matrix and print dimensions
    Eigen::MatrixXd m_noise = sensorNoise(m, g_mean, std_dev, model);
    std::cout << "m_noise: " << m_noise.rows() << " x " << m_noise.cols() << std::endl;

    // apply noise to vector and print dimensions
    Eigen::MatrixXd v_noise = sensorNoise(v, g_mean, std_dev, model);
    std::cout << "v_noise: " << v_noise.rows() << " x " << v_noise.cols() << std::endl;

    return 0;
}