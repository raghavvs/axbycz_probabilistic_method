#include <gtest/gtest.h>
#include <mvg.h>

TEST(mvgTest, testOutputDimensions) {
  Eigen::VectorXd mu(3);
  mu << 1.0, 2.0, 3.0;
  Eigen::MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 1.0;
  int N = 1000;
  Eigen::MatrixXd y, R;
  mvg(mu, Sigma, N, y, R);
  ASSERT_EQ(y.rows(), 3);
  ASSERT_EQ(y.cols(), N);
}

TEST(mvgTest, testSymmetricSigma) {
  Eigen::VectorXd mu(3);
  mu << 1.0, 2.0, 3.0;
  Eigen::MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.5, 0.4, 1.0, 0.5, 0.4, 0.5, 1.0;
  int N = 1000;
  Eigen::MatrixXd y, R;
  ASSERT_DEATH(mvg(mu, Sigma, N, y, R), ".*Sigma must be symmetric\\..*");
}

TEST(mvgTest, testPositiveDefiniteSigma) {
  Eigen::VectorXd mu(3);
  mu << 1.0, 2.0, 3.0;
  Eigen::MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.5, 0.5, 1.0, -0.5, 0.5, -0.5, 1.0;
  int N = 1000;
  Eigen::MatrixXd y, R;
  ASSERT_DEATH(mvg(mu, Sigma, N, y, R), ".*Sigma must be positive definite\\..*");
}

TEST(mvgTest, testZeroSamples) {
  Eigen::VectorXd mu(3);
  mu << 1.0, 2.0, 3.0;
  Eigen::MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 1.0;
  int N = 0;
  Eigen::MatrixXd y, R;
  ASSERT_DEATH(mvg(mu, Sigma, N, y, R), ".*A positive integer number of samples must be requested\\..*");
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}