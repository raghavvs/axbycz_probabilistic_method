#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <skewLog.h>

TEST(SkewLogTest, HandlesIdentityMatrix) {
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d expected = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d result = skewLog(R);
    EXPECT_EQ(result, expected);
}

TEST(SkewLogTest, HandlesArbitraryMatrix) {
    Eigen::Matrix3d R;
    R << 1, 0, 0,
         0, 0.707, -0.707,
         0, 0.707, 0.707;
    Eigen::Matrix3d expected;
    expected << 0, 0, 0,
                0, 0, -1.0,
                0, 1.0, 0;
    Eigen::Matrix3d result = skewLog(R);
    EXPECT_EQ(result, expected);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}