#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>

#include "rotError.h" 

TEST(RotationErrorTest, SameRotationMatrixReturnsZeroError) {
  Eigen::Matrix4d X1, X2;
  X1.setIdentity();
  X2.setIdentity();

  EXPECT_EQ(rotError(X1, X2), 0);
}

TEST(RotationErrorTest, AntiParallelRotationMatricesReturnsPiError) {
  Eigen::Matrix4d X1, X2;
  X1.setIdentity();
  X2 << -1, 0, 0, 0,
         0, -1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;

  EXPECT_EQ(rotError(X1, X2), M_PI);
}

TEST(RotationErrorTest, RandomRotationMatricesReturnsCorrectError) {
  Eigen::Matrix4d X1, X2;
  X1 << -0.684,-0.732,0.0014,-1.54,
        0.718,-0.671,0.1835,3.152,
      -0.1231,-0.124,-0.984,-0.261,
             0,     0,     0,     1;
  X2 << -0.977,0.2062,0.0441,-1.532,
      -0.1976,-0.9473,-0.2513,-3.191,
        0.079,-0.2432,0.9668,-0.261,
             0,     0,     0,     1;

  double expected_error = 0.2826; // pre-computed error

  EXPECT_NEAR(rotError(X1, X2), expected_error, 1e-4); // using a tolerance of 1e-4
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}