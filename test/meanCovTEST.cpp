#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include <../util/expm.h>
#include <../util/meanCov.h>

TEST(MeanCovTest, CalculatesCorrectMeanAndCovariance)
{
    // Generate random rotation matrices and translation vectors
    int N = 10;
    Eigen::MatrixXd X(4, 4 * N);
    for (int i = 0; i < N; i++)
    {
        Eigen::Matrix3d R = Eigen::AngleAxisd(Eigen::Vector3d::Random()).toRotationMatrix();
        Eigen::Vector3d t = Eigen::Vector3d::Random();
        X.block<3, 3>(0, 4 * i) = R;
        X.block<3, 1>(0, 4 * i + 3) = t;
        X.block<1, 4>(3, 4 * i) << 0, 0, 0, 1;
    }

    // Compute mean and covariance
    Eigen::VectorXd Mean;
    Eigen::MatrixXd Cov;
    std::tie(Mean, Cov) = meanCov(X, Mean, Cov);

    // Check that the calculated mean is close to the true mean (should be within 1e-6)
    Eigen::Matrix4d true_mean = Eigen::Matrix4d::Identity();
    for (int i = 0; i < N; i++)
    {
        Eigen::Matrix4d Xi = Eigen::Map<const Eigen::Matrix4d>(X.data() + i * 16, 4, 4);
        true_mean *= Xi;
    }
    true_mean = true_mean.pow(1.0 / N);
    double error = (logm(true_mean) - logm(Mean)).norm();
    EXPECT_LT(error, 1e-6);

    // Check that the calculated covariance is symmetric and positive definite (should be within 1e-6)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Cov);
    EXPECT_TRUE(es.eigenvalues().minCoeff() > 0);
    EXPECT_LT((Cov - Cov.transpose()).norm(), 1e-6);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}