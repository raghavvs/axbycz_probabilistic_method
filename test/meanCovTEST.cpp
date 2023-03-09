#include <gtest/gtest.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <vector>

#include <meanCov.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::Vector3d;

// Test case for meanCov function
/* TEST(MeanCovTest, MeanAndCovariance)
{
    // Generate test data
    const int N = 100;
    MatrixXd X = MatrixXd::Random(4, 4*N);

    // Calculate mean and covariance using the function being tested
    MatrixXd Mean, Cov;
    meanCov(X, Mean, Cov);

    // Calculate mean and covariance using Eigen built-in functions for comparison
    MatrixXd X_mean = MatrixXd::Zero(4, 4);
    for (int i = 0; i < N; i++)
    {
        X_mean += X.block<4, 4>(0, 4*i);
    }
    X_mean /= N;
    MatrixXd X_centered = X.colwise() - X_mean;
    MatrixXd X_cov = (X_centered * X_centered.transpose()) / N;

    // Check that the calculated mean is close to the true mean
    EXPECT_TRUE(Mean.isApprox(X_mean, 1e-5));

    // Check that the calculated covariance is close to the true covariance
    EXPECT_TRUE(Cov.isApprox(X_cov, 1e-5));
} */

TEST(MeanCovTest, TestMeanCov) {
    Eigen::MatrixXd X(4, 16);
    X << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

    Eigen::MatrixXd Mean, Cov;
    meanCov(X, Mean, Cov);

    // Add test assertions here
    ASSERT_EQ(Mean.rows(), 4);
    ASSERT_EQ(Mean.cols(), 4);
    ASSERT_EQ(Cov.rows(), 6);
    ASSERT_EQ(Cov.cols(), 6);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}