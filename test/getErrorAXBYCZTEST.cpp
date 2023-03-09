#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>
#include <getErrorAXBYCZ.h>

/* TEST(getErrorAXBYCZTest, testXYZError) {
  Eigen::Matrix4d X_f = Eigen::Matrix4d::Random();
  Eigen::Matrix4d Y_f = Eigen::Matrix4d::Random();
  Eigen::Matrix4d Z_f = Eigen::Matrix4d::Random();
  Eigen::Matrix4d XActual = Eigen::Matrix4d::Random();
  Eigen::Matrix4d YActual = Eigen::Matrix4d::Random();
  Eigen::Matrix4d ZActual = Eigen::Matrix4d::Random();

  Eigen::Vector3d expectedError;
  expectedError << 0.1, 0.2, 0.3;

  Eigen::Vector3d actualError = getErrorAXBYCZ(X_f, Y_f, Z_f, XActual, YActual, ZActual);
  
  ASSERT_NEAR(expectedError(0), actualError(0), 1e-6);
  ASSERT_NEAR(expectedError(1), actualError(1), 1e-6);
  ASSERT_NEAR(expectedError(2), actualError(2), 1e-6);
}

TEST(getErrorAXBYCZTest, testTranslationError) {
  Eigen::Matrix4d X_f = Eigen::Matrix4d::Random();
  Eigen::Matrix4d Y_f = Eigen::Matrix4d::Random();
  Eigen::Matrix4d Z_f = Eigen::Matrix4d::Random();
  Eigen::Matrix4d XActual = Eigen::Matrix4d::Random();
  Eigen::Matrix4d YActual = Eigen::Matrix4d::Random();
  Eigen::Matrix4d ZActual = Eigen::Matrix4d::Random();

  Eigen::Vector3d expectedError;
  expectedError << 0.01, 0.02, 0.03;

  Eigen::Vector3d actualError = getErrorAXBYCZ(X_f, Y_f, Z_f, XActual, YActual, ZActual);
  
  ASSERT_NEAR(expectedError(3), actualError(3), 1e-6);
  ASSERT_NEAR(expectedError(4), actualError(4), 1e-6);
  ASSERT_NEAR(expectedError(5), actualError(5), 1e-6);
} */

TEST(getErrorAXBYCZTest, Test1) {
  // Initialize random matrices
  Eigen::Matrix4d X_f = Eigen::Matrix4d::Random();
  Eigen::Matrix4d Y_f = Eigen::Matrix4d::Random();
  Eigen::Matrix4d Z_f = Eigen::Matrix4d::Random();
  Eigen::Matrix4d XActual = Eigen::Matrix4d::Random();
  Eigen::Matrix4d YActual = Eigen::Matrix4d::Random();
  Eigen::Matrix4d ZActual = Eigen::Matrix4d::Random();

  // Call function to get error
  Eigen::Vector3d result = getErrorAXBYCZ(X_f, Y_f, Z_f, XActual, YActual, ZActual);

  // Check that the result has the expected dimensions
  ASSERT_EQ(result.rows(), 6);
  ASSERT_EQ(result.cols(), 1);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}