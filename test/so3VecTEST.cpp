#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>

Eigen::Matrix3d so3_vec(const Eigen::Vector3d& X) {
  Eigen::Matrix3d g;
  if (X.size() == 3) {
    g << 0, -X(2), X(1),
         X(2), 0, -X(0),
         -X(1), X(0), 0;
  } else { 
    g << 0, -X(2), X(1),
         X(2), 0, -X(0),
         -X(1), X(0), 0;
  }
  return g;
}

TEST(So3VecTest, HandlesZeroVector) {
  Eigen::Vector3d v(0, 0, 0);
  Eigen::Matrix3d expected;
  expected << 0, 0, 0,
              0, 0, 0,
              0, 0, 0;
  Eigen::Matrix3d result = so3_vec(v);
  EXPECT_EQ(result, expected);
}

TEST(So3VecTest, HandlesNonZeroVector) {
  Eigen::Vector3d v(1, 2, 3);
  Eigen::Matrix3d expected;
  expected << 0, -3, 2,
              3, 0, -1,
              -2, 1, 0;
  Eigen::Matrix3d result = so3_vec(v);
  EXPECT_EQ(result, expected);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}