#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>
#include <expm.h>

TEST(ExpmTest, IdentityMatrix) {
    Eigen::Matrix4d A = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d result = expm(A);
    Eigen::Matrix4d expected = Eigen::Matrix4d::Identity();
    ASSERT_TRUE(result.isApprox(expected));
}

TEST(ExpmTest, ZeroMatrix) {
    Eigen::Matrix4d A = Eigen::Matrix4d::Zero();
    Eigen::Matrix4d result = expm(A);
    Eigen::Matrix4d expected = Eigen::Matrix4d::Identity();
    ASSERT_TRUE(result.isApprox(expected));
}

TEST(ExpmTest, DiagonalMatrix) {
    Eigen::Matrix4d A;
    A << 1.0, 0.0, 0.0, 0.0,
         0.0, 2.0, 0.0, 0.0,
         0.0, 0.0, 3.0, 0.0,
         0.0, 0.0, 0.0, 4.0;
    Eigen::Matrix4d result = expm(A);
    Eigen::Matrix4d expected;
    expected << std::exp(1.0), 0.0, 0.0, 0.0,
                0.0, std::exp(2.0), 0.0, 0.0,
                0.0, 0.0, std::exp(3.0), 0.0,
                0.0, 0.0, 0.0, std::exp(4.0);
    ASSERT_TRUE(result.isApprox(expected));
}