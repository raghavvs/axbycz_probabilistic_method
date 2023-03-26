/*
The code defines a function named fKine that computes the forward kinematics 
of a robotic arm using the Denavit-Hartenberg (DH) convention. The function 
takes in a 6-dimensional vector of joint angles q as input and returns a 4x4 
homogeneous transformation matrix T that describes the pose of the end-effector 
with respect to the robot's base frame. The DH parameters of the robot are 
defined as constants at the beginning of the function, and the function 
iteratively computes the forward kinematics by multiplying the homogeneous 
transformation matrices Ti for each joint together. Finally, the resulting 
homogeneous transformation matrix T is returned.
*/

#ifndef FKINE_H
#define FKINE_H

#include <iostream>
#include <eigen3/Eigen/Dense>

// Define a function to compute the forward kinematics of the robot
Eigen::Matrix4d fKine(const Eigen::VectorXd& q)
{
    // Define the DH parameters of the ABB IRB120 robot (needs to be updated)
    const double d1 = 0.290;
    const double a2 = 0.270;
    const double a3 = 0.070;
    const double d4 = 0.302;
    const double d5 = 0.072;
    const double d6 = 0.074;

    // Define the joint angles (theta)
    Eigen::VectorXd theta(6);
    theta << q(0), q(1), q(2), q(3), q(4), q(5);

    // Define the DH parameters (d, a, alpha)
    Eigen::Matrix<double, 6, 1> d;
    Eigen::Matrix<double, 6, 1> a;
    Eigen::Matrix<double, 6, 1> alpha;
    d << d1, 0, 0, d4, d5, d6;
    a << 0, a2, a3, 0, 0, 0;
    alpha << -M_PI_2, 0, -M_PI_2, M_PI_2, -M_PI_2, 0;

    // Compute the forward kinematics
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    for (int i = 0; i < 6; i++) {
        Eigen::Matrix4d Ti;
            Ti << cos(theta[i]), -sin(theta[i])*cos(alpha[i]),  sin(theta[i])*sin(alpha[i]), a[i]*cos(theta[i]),
              sin(theta[i]),  cos(theta[i])*cos(alpha[i]), -cos(theta[i])*sin(alpha[i]), a[i]*sin(theta[i]),
                       0,             sin(alpha[i]),             cos(alpha[i]),            d[i],
                       0,                       0,                       0,            1;
        T *= Ti;
    }
    return T;
}

#endif