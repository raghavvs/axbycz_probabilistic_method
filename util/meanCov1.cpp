#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

// Define a function that converts a 3x3 skew-symmetric matrix into a 3x1 vector
Eigen::Vector3d vex(Eigen::Matrix3d M) {
    Eigen::Vector3d v;
    v << M(2,1), M(0,2), M(1,0);
    return v;
}

void meanCov(Eigen::MatrixXd* X, int N, Eigen::Matrix4d& Mean, Eigen::MatrixXd& Cov) {
    Mean = Eigen::Matrix4d::Identity();
    Cov = Eigen::MatrixXd::Zero(6,6);
    Eigen::Matrix4d sum_se = Eigen::Matrix4d::Zero();
    for (int i = 0; i < N; i++) {
        sum_se += X[i].log();
    }
    Mean = (sum_se/N).exp();

    // Iterative process to calculate the true Mean
    Eigen::Matrix4d diff_se = Eigen::Matrix4d::Ones();
    int max_num = 100;
    double tol = 1e-5;
    int count = 1;
    while (diff_se.norm() >= tol && count <= max_num) {
        diff_se = Eigen::Matrix4d::Zero();
        for (int i = 0; i < N; i++) {
            diff_se += (Mean.inverse()*X[i]).log();
        }
        Mean *= (diff_se/N).exp();
        count++;
    }

    // Covariance
    for (int i = 0; i < N; i++) {
        diff_se = (Mean.inverse()*X[i]).log();
        Eigen::VectorXd diff_vex(6);
        diff_vex << vex(diff_se.block<3,3>(0,0)), diff_se.block<3,1>(0,3);
        Cov += diff_vex * diff_vex.transpose();
    }
    Cov /= N;
}

int main()
{
    int N = 5;
    Eigen::MatrixXd* X = new Eigen::MatrixXd[N]; // Dynamically allocate memory for the array
    for (int i = 0; i < N; i++) {
        X[i] = Eigen::MatrixXd(4,4); // Initialize each element as a 4x4 matrix
        X[i] << rand()%10+1 , rand()%10+1 , rand()%10+1 , rand()%10+1 , // Fill each element with random numbers from 1 to 10
                rand()%10+1 , rand()%10+1 , rand()%10+1 , rand()%10+1 ,
                rand()%10+1 , rand()%10+1 , rand()%10+1 , rand()%10+1 ,
                rand()%10+1 , rand()%10+1 , rand()%10+1 , rand()%10+1 ;
    }

    // Print the input array
    std::cout << "The input array is: " << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << "X[" << i << "] =" << std::endl;
        std::cout << X[i] << std::endl;
        std::cout << std::endl;
    }

    // Declare variables to store the output mean and covariance
    Eigen::Matrix4d Mean;
    Eigen::MatrixXd Cov;

    // Call the meanCov function with the input and output arguments
    meanCov(X,N,Mean,Cov);

    // Print the output mean and covariance
    std::cout << "The output mean is: " << std::endl;
    std::cout << Mean << std::endl;

    std::cout <<"The output covariance is:"<<std::endl;

    std::cout<<Cov<<std::endl;

    // Free the memory allocated for the input array

    delete[] X;

    return 0;
}
