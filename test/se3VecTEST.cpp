#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>
#include "se3Vec.h" // header file containing the se3Vec function

// Test fixture class for se3Vec function tests
class Se3VecTest : public testing::Test {
protected:
  Eigen::Matrix<double, 4, 4> X1;
  Eigen::Matrix<double, 3, 3> X2;
  Eigen::Matrix<double, 4, 6> X3;

  Se3VecTest() {
    X1 << 1, 2, 3, 4,
          5, 6, 7, 8,
          9, 10, 11, 12,
          13, 14, 15, 16;

    X2 << 1, 2, 3,
          4, 5, 6,
          7, 8, 9;

    X3 << 1, 2, 3, 4, 5, 6,
          7, 8, 9, 10, 11, 12,
          13, 14, 15, 16, 17, 18,
          19, 20, 21, 22, 23, 24;
  }
};

TEST_F(Se3VecTest, Returns6x1MatrixWhenInputIs4x4Matrix) {
  Eigen::Matrix<double, 6, 1> expected_g;
  expected_g << -X1(1, 2), X1(0, 2), -X1(0, 1), X1(0, 3), X1(1, 3), X1(2, 3);

  Eigen::Matrix<double, 6, 1> actual_g = se3Vec(X1);

  ASSERT_EQ(expected_g.rows(), actual_g.rows());
  ASSERT_EQ(expected_g.cols(), actual_g.cols());
}

TEST_F(Se3VecTest, Returns4x4MatrixWhenInputIsNot4x4Matrix) {
  Eigen::Matrix<double, 4, 4> expected_g;
  expected_g << 0, -X2(2, 1), X2(1, 0), 0,
                X2(2, 1), 0, -X2(0, 2), 0,
                -X2(1, 0), X2(0, 2), 0, 0,
                0, 0, 0, 0;

  Eigen::Matrix<double, 4, 4> actual_g = se3Vec(X2);

  ASSERT_EQ(expected_g.rows(), actual_g.rows());
  ASSERT_EQ(expected_g.cols(), actual_g.cols());
}

TEST_F(Se3VecTest, ReturnsCorrectDimensionsForAnotherNon4x4Matrix) {
    
    Eigen::Matrix<double, 4, 4> expected_g;
    expected_g << 0, -X3(2, 1), X3(1, 0), X3(0, 3),
                X3(2, 1), 0,-X3(0, 1), X3(1, 3),
                -X3(1, 0), X3(0, 1), 0, X3(2, 3),
                0, 0, 0, 1;

    // Define the twist vector.
    Eigen::Matrix<double, 6, 1> twist_vec;
    twist_vec << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;

    // Compute the expected result.
    Eigen::Matrix<double, 4, 4> expected_g_twist = expected_g * se3Vec(twist_vec).exp();

    // Compute the actual result.
    Eigen::Matrix<double, 4, 4> actual_g_twist = g(twist_vec);

    // Verify the result.
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ASSERT_NEAR(expected_g_twist(i, j), actual_g_twist(i, j), 1e-9);
         }
     }
}