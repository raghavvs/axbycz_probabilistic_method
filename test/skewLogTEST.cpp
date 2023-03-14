#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <skewLog.h>

Eigen::Matrix3d skewLog(Eigen::Matrix3d R);

TEST(SkewLogTest, HandlesIdentityMatrix) {
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d result = skewLog(I);
  Eigen::Matrix3d expected = Eigen::Matrix3d::Zero();
  ASSERT_TRUE(result.isApprox(expected));
}

TEST(SkewLogTest, HandlesArbitraryMatrix) {
  Eigen::Matrix3d R;
  R << 0, 0, 1,
       0, -1, 0,
       1, 0, 0;
  Eigen::Matrix3d result = skewLog(R);
  Eigen::Matrix3d expected;
  expected << 0, 0, 0,
              0, 0, -std::atan(1.0),
              0, std::atan(1.0), 0;
  ASSERT_TRUE(result.isApprox(expected, 1e-6));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}