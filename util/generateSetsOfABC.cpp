/*
This generates sets of three vectors, each consisting of a 
constant and two free variables, namely A1, B1, C1, A2, B2, C2, 
A3, B3, C3. The function takes as input the number of sets to 
generate (Num), an option for the probability density function 
(optPDF), a mean vector (Mean), a covariance matrix (Cov), and 
actual X, Y, and Z vectors (XActual, YActual, ZActual). 

The function uses the "generateABC"to generate the sets of ABC values. 
The first set of vectors has a constant value for A1, the second 
set has a constant value for C2, and the third set has a constant 
value for B3. The free variables in each set are randomly generated 
using the specified options and input parameters.
*/

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <generateABC.h>

void generateSetsOfABC(int Num, int optPDF, Eigen::Vector3d Mean, Eigen::Matrix3d Cov, 
                        Eigen::VectorXd XActual, Eigen::VectorXd YActual, Eigen::VectorXd ZActual, 
                        Eigen::Matrix4d& A1, Eigen::Matrix4d& B1, Eigen::Matrix4d& C1, Eigen::Matrix4d& A2, 
                        Eigen::Matrix4d& B2, Eigen::Matrix4d& C2, Eigen::Matrix4d& A3, Eigen::Matrix4d& B3, Eigen::Matrix4d& C3) {

    // Generate constant A1, free B1 and C1
    int opt = 1;
    generateABC(Num, opt, optPDF, Mean, Cov, XActual, YActual, ZActual, A1, B1, C1);

    // Generate constant C2, free A2 and B2
    opt = 3;
    generateABC(Num, opt, optPDF, Mean, Cov, XActual, YActual, ZActual, A2, B2, C2);

    // Generate constant B3, free A3 and C3
    opt = 2;
    generateABC(Num, opt, optPDF, Mean, Cov, XActual, YActual, ZActual, A3, B3, C3);
}