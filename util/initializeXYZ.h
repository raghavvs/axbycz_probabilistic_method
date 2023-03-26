/*
DESCRIPTION:

The program InitializeXYZ generates a set of 3 matrices X, Y, and Z, 
based on the input option opt. If opt is 1, it generates 3 random 
matrices using the randn function and transforms them into 4x4 
homogeneous transformation matrices using the expm function. If 
opt is 2, it sets the matrices X, Y, and Z to pre-defined values. 
If opt is 3, it sets the matrices X, Y, and Z to pre-defined rotation 
matrices. If opt is not 1, 2, or 3, it outputs an error message. 
The matrices X, Y, and Z represent the ground truth transformations 
for a 3D pose estimation problem.
*/

#ifndef INITIALIZEXYZ_H
#define INITIALIZEXYZ_H

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <se3Vec.h>
#include <expm.h>

void initializeXYZ(int opt,
                   Eigen::Matrix4d& X,
                   Eigen::Matrix4d& Y,
                   Eigen::Matrix4d& Z)
{
    if (opt == 1)
    {
        Eigen::VectorXd x = Eigen::VectorXd::Random(6);
        x.normalize();
        X = Eigen::Matrix4d::Identity() * expm(se3Vec(x));
        
        Eigen::VectorXd y = Eigen::VectorXd::Random(6);
        y.normalize();
        Y = Eigen::Matrix4d::Identity() * expm(se3Vec(y));
        
        Eigen::VectorXd z = Eigen::VectorXd::Random(6);
        z.normalize();
        Z = Eigen::Matrix4d::Identity() * expm(se3Vec(z));
    }
    else if (opt == 2)
    {
        X << -0.9765,  0.0636, -0.2059,  0.0215,
              -0.0947, -0.9849,  0.1447, -0.0029,
              -0.1936,  0.1608,  0.9678, -0.0597,
                   0.0,       0.0,       0.0,    1.0;
        
        Y << -0.99908, -0.03266,  0.02786,  164.226/1000,
               0.02737,  0.01553,  0.99950,  301.638/1000,
              -0.03308,  0.99935, -0.01462, -962.841/1000,
                    0.0,       0.0,       0.0,        1.0;
        
        Z <<  0.70063, -0.40451,  0.58779, 0.006,
               0.69084,  0.17849, -0.70063, 0.030,
               0.17849,  0.89695,  0.40451, 0.921,
                    0.0,       0.0,       0.0, 1.0;
    }
    else if (opt == 3)
    {
        X.block<3, 3>(0, 0) = Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()).toRotationMatrix() *
                              Eigen::AngleAxisd(M_PI/3, Eigen::Vector3d::UnitX()).toRotationMatrix();
        
        Y.block<3, 3>(0, 0) = Eigen::AngleAxisd(M_PI/6, Eigen::Vector3d::UnitY()).toRotationMatrix() *
                              Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        
        Z.block<3, 3>(0, 0) = Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()).toRotationMatrix() *
                              Eigen::AngleAxisd(M_PI/3, Eigen::Vector3d::UnitX()).toRotationMatrix() *
                              Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitY()).toRotationMatrix();
    }
    else
    {
        std::cout << "The given opt = " << opt << " is not an option." << std::endl;
        return;
    }
}

#endif