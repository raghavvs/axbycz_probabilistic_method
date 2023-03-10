#include <gtest/gtest.h>
#include "metric.h"

TEST(MetricTest, Test1) {
    std::vector<Eigen::Matrix3d> A(2);
    std::vector<Eigen::Matrix3d> B(2);
    std::vector<Eigen::Matrix3d> C(2);

    A[0] << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    A[1] << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    B[0] << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    B[1] << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    C[0] << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    C[1] << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    Eigen::Matrix3d X, Y, Z;
    X << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    Y << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    Z << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    double result = metric(A, B, C, X, Y, Z);

    ASSERT_NEAR(result, 20.7521, 1e-4);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}