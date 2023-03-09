/*
DESCRIPTION:

This program defines several functions for calculating the mean 
and covariance of a set of 4x4 matrices. The log function takes 
a 4x4 matrix and calculates its matrix logarithm, while the vex 
function takes a 3x3 matrix and returns its vector of exteriorization. 
The meanCov function takes an array of 3x3 matrices and its size 
N and calculates the mean and covariance of the logarithms of those 
matrices. It does this by first taking the average of the logarithms 
using the expm function, then iteratively refining this average until 
convergence using the log and vex functions. Finally, it calculates 
the covariance by taking the vector of differences between each logarithm 
and the mean logarithm and computing their outer product.
*/

#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <vector>

Eigen::Vector3d vex(const Eigen::Matrix3d& S)
{
    Eigen::Vector3d w;
    w << S(2, 1), S(0, 2), S(1, 0);
    return w;
}

void meanCov(Eigen::MatrixXd X, Eigen::MatrixXd& Mean, Eigen::MatrixXd& Cov)
{
    int N = X.cols() / 4;
    Mean = Eigen::MatrixXd::Identity(4, 4);
    Cov = Eigen::MatrixXd::Zero(6, 6);

    // Initial approximation of Mean
    Eigen::MatrixXd sum_se = Eigen::MatrixXd::Zero(4, 4);
    for (int i = 0; i < N; i++)
    {
        sum_se += X.block<4, 4>(0, 4*i).log();
    }
    Mean = (sum_se / N).exp();

    // Iterative process to calculate the true Mean
    Eigen::MatrixXd diff_se = Eigen::MatrixXd::Ones(4, 4);
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num)
    {
        diff_se.setZero();
        for (int i = 0; i < N; i++)
        {
            diff_se += (Mean.inverse() * X.block<4, 4>(0, 4*i)).log();
        }
        Mean = Mean * diff_se.exp() / N;
        count++;
    }

    // Covariance
    for (int i = 0; i < N; i++)
    {
        Eigen::MatrixXd diff_se = (Mean.inverse() * X.block<4, 4>(0, 4*i)).log();
        Eigen::VectorXd diff_vex(6);
        diff_vex << vex(diff_se.block<3, 3>(0, 0)), diff_se.block<3, 1>(0, 3);
        Cov += diff_vex * diff_vex.transpose();
    }
    Cov /= N;
}

/* int main()
{
    std::cout << "Build successful" << std::endl;
    return 0;
} */

/*
vex() is a function that takes a skew-symmetric matrix 
and returns its corresponding 3D vector.
*/

//TEST Case

/* int main()
{
    // Define some 4x4 matrices
    Eigen::MatrixXd X1;
    X1 << 1, 0, 0, 1,
          0, 1, 0, 2,
          0, 0, 1, 3,
          0, 0, 0, 1;
          
    Eigen::MatrixXd X2;
    X2 << -1, 0, 0, 4,
          0, -1, 0, 5,
          0, 0, -1, 6,
          0, 0, 0, 1;
          
    Eigen::MatrixXd X3;
    X3 << 0, 1, 0, 7,
          -1, 0, 0, 8,
          0, 0, 1, 9,
          0, 0, 0, 1;

    // Put the matrices in a vector
    Eigen::MatrixXd X(16, 3);
    Eigen::Map<Eigen::MatrixXd>(X.data(), 4, 4) = X1;
    Eigen::Map<Eigen::MatrixXd>(X.data() + 16, 4, 4) = X2;
    Eigen::Map<Eigen::MatrixXd>(X.data() + 32, 4, 4) = X3;

    // Compute the mean and covariance
    Eigen::VectorXd Mean(4);
    Eigen::MatrixXd Cov(6, 6);
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> result = meanCov(X, Mean, Cov);

    // Print the results
    std::cout << "Mean:\n" << result.first << "\n";
    std::cout << "Covariance:\n" << result.second << "\n";

    return 0;
} */