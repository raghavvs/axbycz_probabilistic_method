#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include "se3Vec.h"

TEST(se3VecTest, MatrixConversion) {
    Eigen::Matrix4d inputMatrix;
    inputMatrix << 1, 0, 0, 1,
                   0, 1, 0, 2,
                   0, 0, 1, 3,
                   0, 0, 0, 1;

    Eigen::Matrix<double, 6, 1> expectedOutput;
    expectedOutput << 0, 0, 0, 1, 2, 3;

    Eigen::Matrix<double, 6, 1> actualOutput = se3Vec(inputMatrix);

    ASSERT_EQ(actualOutput, expectedOutput);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}