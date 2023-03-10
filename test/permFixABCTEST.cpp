#include <gtest/gtest.h>
#include <permFixABC.h>

TEST(PermFixABCTest, OutputMatricesHaveCorrectDimensions) {
  Eigen::MatrixXd M = Eigen::MatrixXd::Random(3, 4);
  Eigen::MatrixXd N = Eigen::MatrixXd::Random(5, 6);
  Eigen::MatrixXd P = Eigen::MatrixXd::Random(7, 8);
  double r = 0.5;
  Eigen::MatrixXd M_perm, N_perm, P_perm;

  permFixABC(M, N, P, r, M_perm, N_perm, P_perm);

  EXPECT_EQ(M_perm.rows(), M.rows());
  EXPECT_EQ(M_perm.cols(), M.cols() * N.cols());
  EXPECT_EQ(N_perm.rows(), N.rows());
  EXPECT_EQ(N_perm.cols(), N.cols());
  EXPECT_EQ(P_perm.rows(), P.rows());
  EXPECT_EQ(P_perm.cols(), P.cols());
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}