#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <se3Vec.h>
#include <expm.h>

void InitializeXYZ(int opt, Eigen::Matrix4d& X, Eigen::Matrix4d& Y, Eigen::Matrix4d& Z);

// Test fixture for InitializeXYZ function
class InitializeXYZTest : public ::testing::Test {
 protected:
  void SetUp() override {
    opt_ = 1;
  }

  // Variables to be used in the tests
  int opt_;
  Eigen::Matrix4d X_, Y_, Z_;
};

// Test that InitializeXYZ function initializes matrices correctly for opt = 1
TEST_F(InitializeXYZTest, TestOpt1) {
  InitializeXYZ(opt_, X_, Y_, Z_);
  
  // Check if the matrices have the correct dimensions
  EXPECT_EQ(X_.rows(), 4);
  EXPECT_EQ(X_.cols(), 4);
  EXPECT_EQ(Y_.rows(), 4);
  EXPECT_EQ(Y_.cols(), 4);
  EXPECT_EQ(Z_.rows(), 4);
  EXPECT_EQ(Z_.cols(), 4);
  
  // Check if the matrices are invertible
  EXPECT_NEAR(X_.determinant(), 1.0, 1e-10);
  EXPECT_NEAR(Y_.determinant(), 1.0, 1e-10);
  EXPECT_NEAR(Z_.determinant(), 1.0, 1e-10);
}

// Test that InitializeXYZ function initializes matrices correctly for opt = 2
TEST_F(InitializeXYZTest, TestOpt2) {
  opt_ = 2;
  InitializeXYZ(opt_, X_, Y_, Z_);
  
  // Check if the matrices have the correct values
  Eigen::Matrix4d X_expected, Y_expected, Z_expected;
  X_expected << -0.9765,  0.0636, -0.2059,  0.0215,
                -0.0947, -0.9849,  0.1447, -0.0029,
                -0.1936,  0.1608,  0.9678, -0.0597,
                     0.0,       0.0,       0.0,    1.0;
  Y_expected << -0.99908, -0.03266,  0.02786,  164.226/1000,
                 0.02737,  0.01553,  0.99950,  301.638/1000,
                -0.03308,  0.99935, -0.01462, -962.841/1000,
                      0.0,       0.0,       0.0,        1.0;
  Z_expected <<  0.70063, -0.40451,  0.58779, 0.006,
                 0.69084,  0.17849, -0.70063, 0.030,
                 0.17849,  0.89695,  0.40451, 0.921,
                      0.0,       0.0,       0.0, 1.0;
  EXPECT_TRUE(X_.isApprox(X_expected));
  EXPECT_TRUE(Y_.isApprox(Y_expected));
  EXPECT_TRUE(Z_.isApprox(Z_expected));
}

// Test that InitializeXYZ function initializes matrices correctly for opt = 3
TEST_F(InitializeXYZTest, TestOpt3) {
    opt_ = 3;
    InitializeXYZ(opt_, X_, Y_, Z_);

    Eigen::Matrix4d X_expected, Y_expected, Z_expected;
    X_expected.setIdentity();
    Y_expected.setIdentity();
    Z_expected.setIdentity();
    X_expected.block<3, 3>(0, 0) = Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()).toRotationMatrix() *
    Eigen::AngleAxisd(M_PI/3, Eigen::Vector3d::UnitX()).toRotationMatrix();
    Y_expected.block<3, 3>(0, 0) = Eigen::AngleAxisd(M_PI/6, Eigen::Vector3d::UnitY()).toRotationMatrix() *
    Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    Z_expected.block<3, 3>(0, 0) = Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()).toRotationMatrix() *
    Eigen::AngleAxisd(M_PI/3, Eigen::Vector3d::UnitX()).toRotationMatrix() *
    Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitY()).toRotationMatrix();

    EXPECT_NEAR(X_.determinant(), 1.0, 1e-10);
    EXPECT_NEAR(Y_.determinant(), 1.0, 1e-10);
    EXPECT_NEAR(Z_.determinant(), 1.0, 1e-10);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}