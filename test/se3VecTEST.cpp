#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include "se3Vec.h"

TEST(Se3VecTest, Handles4x4Matrix) {
    Eigen::Matrix4d X;
    X << 1, 0, 0, 2,
         0, 1, 0, 3,
         0, 0, 1, 4,
         0, 0, 0, 1;

    Eigen::Matrix<double, 6, 1> expected;
    expected << 0, 0, 0, 2, 3, 4;

    Eigen::Matrix<double, 6, 1> result = se3Vec(X);

    ASSERT_TRUE(result.isApprox(expected));
}

TEST(Se3VecTest, Handles4x3Matrix) {
    Eigen::Matrix<double, 4, 3> X;
    X << 1, 0, 0,
         0, 1, 0,
         0, 0, 1,
         2, 3, 4;

    Eigen::Matrix<double, 6, 1> expected;
    expected << 0, -4, 3, 2, 4, -1;

    Eigen::Matrix<double, 6, 1> result = se3Vec(X);

    ASSERT_TRUE(result.isApprox(expected));
}



int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}