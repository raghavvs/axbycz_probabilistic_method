#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>
#include "scrambleData.h"

TEST(ScrambleDataTest, ReturnsSameMatrixWhenScrambleRateIsZero) {
    // Create input matrix
    Eigen::MatrixXd M(3, 3);
    M << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    // Call function with zero scramble rate
    Eigen::MatrixXd M_scrambled = scrambleData(M, 0.0);

    // Check that the output is the same as the input
    ASSERT_TRUE(M_scrambled.isApprox(M));
}

TEST(ScrambleDataTest, ReturnsDifferentMatrixWhenScrambleRateIsNonZero) {
    // Create input matrix
    Eigen::MatrixXd M(3, 3);
    M << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    // Call function with non-zero scramble rate
    Eigen::MatrixXd M_scrambled = scrambleData(M, 0.5);

    // Check that the output is different from the input
    ASSERT_FALSE(M_scrambled.isApprox(M));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}