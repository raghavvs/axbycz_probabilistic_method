#include <cmath>
#include <iostream>
#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>

Eigen::MatrixXd expm(const Eigen::MatrixXd& A) 
{
    Eigen::MatrixXd expA = Eigen::MatrixXd::Identity(A.rows(), A.cols());
    double factorial = 1;
    int i = 1;
    Eigen::MatrixXd An = A;
    
    do
    {
        factorial *= i;
        expA += An / factorial;
        An = An * A;
        ++i;
    } while ((An / factorial).lpNorm<Eigen::Infinity>() > 1e-15);
  
    return expA;
}

TEST(expmTest, MatrixExponential) {
    Eigen::MatrixXd A(3, 3);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    Eigen::MatrixXd result = expm(A);

    Eigen::MatrixXd expected_result(3, 3);
    expected_result << 1118906.69941320,    1374815.06293582,   1630724.42645844,
                        2533881.04189899,    3113415.03138058,    3692947.02086217,
                        3948856.38438479,    4852012.99982534,    5755170.61526590;


    ASSERT_TRUE(result.isApprox(expected_result, 1e-4));

    // Return result
    SUCCEED() << "\nResult:\n" << result;
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    std::cout << "Running main() from expmTEST.cpp" << std::endl;
    RUN_ALL_TESTS();

    // Access the result variable here
    Eigen::MatrixXd A(3, 3); 
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    Eigen::MatrixXd result = expm(A);
    std::cout << "\nResult from main():\n" << result;
    return 0;
}