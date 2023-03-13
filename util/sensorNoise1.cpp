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

                g_noise.col(i) = g.col(i) * (MatrixXd((se3Vec(noise_old1))).exp() * MatrixXd((se3Vec(noise_old2))).exp());
            }
            break;

        default:
            std::cerr << "Invalid noise model specified" << std::endl;
            break;
    }

    return g_noise;
}

int main() {
    // Create a random 4x4 matrix
    MatrixXd g = MatrixXd::Random(4, 4);

    // Set the mean and standard deviation for the noise
    MatrixXd gmean = MatrixXd::Zero(4, 4);
    double std = 0.1;

    // Add noise to the matrix using model 1
    MatrixXd g_noise = sensorNoise(g, gmean, std, 1);

    // Print the original and noisy matrices
    std::cout << "Original matrix:" << std::endl << g << std::endl << std::endl;
    std::cout << "Noisy matrix:" << std::endl << g_noise << std::endl;

    return 0;
}