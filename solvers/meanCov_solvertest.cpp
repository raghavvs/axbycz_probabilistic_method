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

#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "meanCov.h"

int main()
{
    int N = 5;
    //Eigen::MatrixXd* X = new Eigen::MatrixXd[N]; // Dynamically allocate memory for the array
    /*for (int i = 0; i < N; i++) {
        X[i] = Eigen::Matrix4d::Random();
    }*/

    Eigen::Matrix4d X[N];
    for (int i = 0; i < N; i++) {
        X[i] = Eigen::Matrix4d::Random();
    }

    // Print the input array
    std::cout << "The input array is: " << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << "X[" << i << "] =" << std::endl;
        std::cout << X[i] << std::endl;
        std::cout << std::endl;
    }

    // Declare variables to store the output mean and covariance
    Eigen::MatrixXd Mean;
    Eigen::MatrixXd Cov;

    // Call the meanCov function with the input and output arguments
    meanCov(X,N,Mean,Cov);

    // Print the output mean and covariance
    std::cout << "The output mean is: " << std::endl;
    std::cout << Mean << std::endl;

    std::cout <<"The output covariance is:"<<std::endl;

    std::cout<<Cov<<std::endl;

    // Free the memory allocated for the input array

    //delete[] X;

    return 0;
}