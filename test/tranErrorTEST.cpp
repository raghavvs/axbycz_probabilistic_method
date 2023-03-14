#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>
#include "tranError.h"

TEST(tranErrorTest, ComputesCorrectTranslationError) {
  Eigen::MatrixXd X1(4, 4);
  X1 << 1, 0, 0, 2,
        0, 1, 0, 3,
        0, 0, 1, 4,
        0, 0, 0, 1;

  Eigen::MatrixXd X2(4, 4);
  X2 << 1, 0, 0, 5,
        0, 1, 0, 6,
        0, 0, 1, 7,
        0, 0, 0, 1;

  double expected = 5.19615;  // expected translation error
  double actual = tranError(X1, X2);

  EXPECT_NEAR(actual, expected, 1e-5);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}