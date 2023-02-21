#include <eigen3/Eigen/Dense>
#include <iostream>
#include <random>

using namespace Eigen;

// Function to convert a 6-vector to a 4x4 matrix
Matrix4d se3_vec(const VectorXd& xi)
{
    Matrix4d mat;
    mat << 0, -xi(5), xi(4), xi(0),
           xi(5), 0, -xi(3), xi(1),
           -xi(4), xi(3), 0, xi(2),
           0, 0, 0, 0;
    return mat.exp();
}

void InitializeXYZ(int opt, Matrix4d& X, Matrix4d& Y, Matrix4d& Z)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1);
    VectorXd x(6), y(6), z(6);

    if (opt == 1)
    {
        // Generate a random X
        for (int i = 0; i < 6; ++i)
            x(i) = d(gen);
        x /= x.norm();
        X = se3_vec(x);

        // Generate a random Y
        for (int i = 0; i < 6; ++i)
            y(i) = d(gen);
        y /= y.norm();
        Y = se3_vec(y);

        // Generate a random Z
        for (int i = 0; i < 6; ++i)
            z(i) = d(gen);
        z /= z.norm();
        Z = se3_vec(z);
    }
    else if (opt == 2)
    {
        X << -0.9765, 0.0636, -0.2059, 0.0215,
             -0.0947, -0.9849, 0.1447, -0.0029,
             -0.1936, 0.1608, 0.9678, -0.0597,
             0, 0, 0, 1;

        Y << -0.99908, -0.03266, 0.02786, 164.226/1000,
             0.02737, 0.01553, 0.99950, 301.638/1000,
             -0.03308, 0.99935, -0.01462, -962.841/1000,
             0, 0, 0, 1;

        Z << 0.70063, -0.40451, 0.58779, 0.006,
             0.69084, 0.17849, -0.70063, 0.030,
             0.17849, 0.89695, 0.40451, 0.921,
             0, 0, 0, 1;
    }
    else if (opt == 3)
    {
        X = AngleAxisd(M_PI/4, Vector3d::UnitZ()) * AngleAxisd(M_PI/3, Vector3d::UnitX());
        Y = AngleAxisd(M_PI/6, Vector3d::UnitY()) * AngleAxisd(M_PI/4, Vector3d::UnitZ());
        Z = AngleAxisd(M_PI/4, Vector3d::UnitY()) * AngleAxisd(M_PI/3, Vector3d::UnitZ()) * AngleAxisd(M_PI)

    }
    else 
    {
        std::cout<<"The given opt = %d is not an option. "<<opt<<std::endl;
        return
    }
}