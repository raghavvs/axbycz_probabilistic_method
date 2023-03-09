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

int main() 
{
    Eigen::MatrixXd X(4, 16);
    X << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

    Eigen::MatrixXd Mean, Cov;
    meanCov(X, Mean, Cov);

    std::cout << "Mean: " << std::endl << Mean << std::endl;
    std::cout << "Covariance: " << std::endl << Cov << std::endl;

     return 0;
}

/* int main() 
{
    MatrixXd X(4, 4 * 3);
    X << 1, 0, 0, 1,
         0, 0, 1, 1,
         0, 1, 0, 1,
         1, 1, 1, 0,
         0, 1, 1, 1,
         1, 0, 1, 1,
         1, 1, 0, 1,
         1, 1, 1, 1,
         0, 1, 1, 1,
         1, 0, 1, 1,
         1, 1, 0, 1,
         1, 1, 1, 1;

    MatrixXd Mean;
    MatrixXd Cov;
    meanCov(X, Mean, Cov);

    std::cout << "Mean: " << std::endl << Mean << std::endl;
    std::cout << "Covariance: " << std::endl << Cov << std::endl;

     return 0;
} */

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
    result = meanCov(X, Mean, Cov);

    // Print the results
    std::cout << "Mean:\n" << result.first << "\n";
    std::cout << "Covariance:\n" << result.second << "\n";

    return 0;
} */

/* int main()
{
    // Generate some random 4x4 matrices
    int N = 10;
    MatrixXd X(4, 4*N);
    for (int i = 0; i < N; i++)
    {
        Matrix4d M = Matrix4d::Random();
        X.block<4, 4>(0, 4*i) = M.exp();
    }

    // Calculate the mean and covariance
    MatrixXd Mean, Cov;
    meanCov(X, Mean, Cov);

    // Print the results
    cout << "Mean: " << endl << Mean << endl;
    cout << "Covariance: " << endl << Cov << endl;

    return 0;
} */

/* int main()
{
    // Define random dataset
    int N = 10; // number of samples
    MatrixXd X(4, 4*N); // data matrix
    for (int i = 0; i < N; i++) {
        Matrix4d T = Matrix4d::Identity();
        Quaterniond q = Quaterniond(AngleAxisd(rand() % 10 / 10.0, Vector3d::UnitX())
                              * AngleAxisd(rand() % 10 / 10.0, Vector3d::UnitY())
                              * AngleAxisd(rand() % 10 / 10.0, Vector3d::UnitZ()));
        Matrix3d R = q.toRotationMatrix();

        T.block<3, 3>(0, 0) = R;
        T(0, 3) = rand() % 10;
        T(1, 3) = rand() % 10;
        T(2, 3) = rand() % 10;
        X.block<4, 4>(0, 4*i) = T;
    }

    // Calculate mean and covariance
    Matrix4d Mean;
    MatrixXd Cov;
    meanCov(X, Mean, Cov);

    // Print results
    cout << "Mean:" << endl << Mean << endl << endl;
    cout << "Covariance:" << endl << Cov << endl;

    return 0;
} */