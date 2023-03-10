#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>
#include "skewExp.h"

TEST(SkewTest, ReturnsCorrectMatrix) {
  Eigen::Vector3d v(1, 2, 3);
  Eigen::Matrix3d expected_m;
  expected_m << 0, -3, 2,
                3, 0, -1,
                -2, 1, 0;
  Eigen::Matrix3d actual_m = skew(v);
  ASSERT_TRUE(actual_m.isApprox(expected_m));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}