#include <gtest/gtest.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <random>
#include <algorithm>
#include <chrono>
#include "scrambleData.h" 

TEST(ScrambleDataTest, OutputMatrixHasSameShape) {
    Eigen::MatrixXd M(3, 3);
    M << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    double s_rate = 1.0;
    Eigen::MatrixXd M_perm = scrambleData(M, s_rate);
    ASSERT_EQ(M_perm.rows(), M.rows());
    ASSERT_EQ(M_perm.cols(), M.cols());
}


TEST(ScrambleDataTest, OutputMatrixIsScrambled) {
    // Test that the output matrix is actually scrambled
    Eigen::MatrixXd inputMatrix(3, 3);
    inputMatrix << 1, 2, 3,
                   4, 5, 6,
                   7, 8, 9;
    Eigen::MatrixXd outputMatrix = scrambleData(inputMatrix, 0.5);
    bool isSame = true;
    for (int i = 0; i < inputMatrix.cols(); i++) {
        if (inputMatrix.col(i) != outputMatrix.col(i)) {
            isSame = false;
            break;
        }
    }
    EXPECT_FALSE(isSame);
}

TEST(ScrambleDataTest, OutputMatrixWithZeroScrambleRate) {
    // Test that the output matrix is the same as the input matrix if scramble rate is 0
    Eigen::MatrixXd inputMatrix(3, 3);
    inputMatrix << 1, 2, 3,
                   4, 5, 6,
                   7, 8, 9;
    Eigen::MatrixXd outputMatrix = scrambleData(inputMatrix, 0.0);
    EXPECT_EQ(inputMatrix, outputMatrix);
}

// Run the test case
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}