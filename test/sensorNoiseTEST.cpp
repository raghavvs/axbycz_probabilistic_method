#include <gtest/gtest.h>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include "se3Vec.h"
#include "so3Vec.h"
#include "sensorNoise.h"

// Include the function being tested
std::vector<Eigen::MatrixXd> sensorNoise(const std::vector<Eigen::MatrixXd> &g, const Eigen::MatrixXd &gmean, const double &std, const int &model);

// Define the test case
TEST(SensorNoiseTest, NoiseModel1) {
    // Set up the test inputs
    Eigen::MatrixXd g1 = Eigen::MatrixXd::Identity(4, 4);
    Eigen::MatrixXd g2 = Eigen::MatrixXd::Identity(4, 4);
    std::vector<Eigen::MatrixXd> g_vec;
    g_vec.push_back(g1);
    g_vec.push_back(g2);
    Eigen::MatrixXd gmean = Eigen::MatrixXd::Zero(6, 1);
    double std = 0.1;
    int model = 1;

    // Compute the expected output
    std::vector<Eigen::MatrixXd> expected_g_noise;
    for (int i = 0; i < g_vec.size(); i++) {
        Eigen::VectorXd temp = Eigen::VectorXd::Random(3);

        Eigen::VectorXd noise_old1 = gmean + Eigen::VectorXd::Zero(6);
        noise_old1.tail(3) = std * Eigen::VectorXd::Random(3);

        Eigen::VectorXd noise_old2 = gmean + Eigen::VectorXd::Zero(6);
        noise_old2.head(3) = std * (temp / temp.norm());
        
        Eigen::MatrixXd g_temp = g_vec[i] * (se3Vec(noise_old1)).exp() * (se3Vec(noise_old2)).exp();
        expected_g_noise.push_back(g_temp);
    }

    // Call the function being tested
    std::vector<Eigen::MatrixXd> actual_g_noise = sensorNoise(g_vec, gmean, std, model);

    // Check that the actual output matches the expected output
    ASSERT_EQ(expected_g_noise.size(), actual_g_noise.size());
    for (int i = 0; i < expected_g_noise.size(); i++) {
        ASSERT_TRUE(expected_g_noise[i].isApprox(actual_g_noise[i], 1e-6));
    }
}

// Run the test case
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}