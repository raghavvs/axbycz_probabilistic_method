#include <gtest/gtest.h>
#include <isEqualf.h>

TEST(isEqualfTest, ReturnsTrueForEqualValues)
{
    double a = 1.0;
    double b = 1.0;

    bool result = isEqualf(a, b);

    EXPECT_TRUE(result);
}

TEST(isEqualfTest, ReturnsFalseForDifferentValues)
{
    double a = 1.0;
    double b = 2.0;

    bool result = isEqualf(a, b);

    EXPECT_FALSE(result);
}

TEST(isEqualfTest, ReturnsFalseForNaNValues)
{
    double a = std::numeric_limits<double>::quiet_NaN();
    double b = 1.0;

    bool result = isEqualf(a, b);

    EXPECT_FALSE(result);
}

TEST(isEqualfTest, ReturnsTrueForValuesWithinThreshold)
{
    double a = 1.0;
    double b = 1.0 + std::numeric_limits<double>::epsilon()*50;

    bool result = isEqualf(a, b);

    EXPECT_TRUE(result);
}

TEST(isEqualfTest, ReturnsFalseForValuesOutsideThreshold)
{
    double a = 1.0;
    double b = 1.0 + std::numeric_limits<double>::epsilon()*150;

    bool result = isEqualf(a, b);

    EXPECT_FALSE(result);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}