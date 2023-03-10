#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>
#include "fKine.h"  


TEST(FKineTest, ComputesCorrectTransform) {
  // Define input joint angles
  Eigen::VectorXd q(6);
  q << 0, 0, 0, 0, 0, 0;

  // Define expected output transform
  Eigen::Matrix4d expected_T;
  expected_T << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;

  // Compute actual output transform using fKine function
  Eigen::Matrix4d actual_T = fKine(q);

  // Check if actual output matches expected output within a tolerance
  EXPECT_TRUE(actual_T.isApprox(expected_T, 1e-6));
}


int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}