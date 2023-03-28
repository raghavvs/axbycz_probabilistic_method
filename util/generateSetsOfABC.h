//
// Created by Raghavendra N S on 3/16/23.
//

/*
DESCRIPTION:

This generates sets of three vectors, each consisting of a
constant and two free variables, namely A1, B1, C1, A2, B2, C2,
A3, B3, C3. The function takes as input the number of sets to
generate (Num), an option for the probability density function
(optPDF), a mean vector (Mean), a covariance matrix (Cov), and
actual X, Y, and Z vectors (XActual, YActual, ZActual).

The function uses the "generateABC" to generate the sets of ABC values.
The first set of vectors has a constant value for A1, the second
set has a constant value for C2, and the third set has a constant
value for B3. The free variables in each set are randomly generated
using the specified options and input parameters.
*/

#ifndef GENERATESETSOFABC_H
#define GENERATESETSOFABC_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <tuple>
#include <vector>
#include "generateABC.h"

struct ABCSets {
    std::vector<Eigen::Matrix4d> A1;
    std::vector<Eigen::Matrix4d> B1;
    std::vector<Eigen::Matrix4d> C1;
    std::vector<Eigen::Matrix4d> A2;
    std::vector<Eigen::Matrix4d> B2;
    std::vector<Eigen::Matrix4d> C2;
    std::vector<Eigen::Matrix4d> A3;
    std::vector<Eigen::Matrix4d> B3;
    std::vector<Eigen::Matrix4d> C3;
};

ABCSets generateSetsOfABC(int Num,
                          int optPDF,
                          const Eigen::VectorXd& Mean,
                          const Eigen::MatrixXd& Cov,
                          const Eigen::Matrix4d& XActual,
                          const Eigen::Matrix4d& YActual,
                          const Eigen::Matrix4d& ZActual) {

    ABCSets result;

    // Generate constant A1, free B1 and C1
    int opt = 1;
    std::tie(result.A1,result.B1,result.C1) = generateABC(Num,opt,optPDF,Mean,Cov,XActual,YActual,ZActual);

    // Generate constant C2, free A2 and B2
    opt = 3;
    std::tie(result.A2,result.B2,result.C2) = generateABC(Num,opt,optPDF,Mean,Cov,XActual,YActual,ZActual);

    // Generate constant B3, free A3 and C3
    opt = 2;
    std::tie(result.A3,result.B3,result.C3) = generateABC(Num,opt,optPDF,Mean,Cov,XActual,YActual,ZActual);

    return result;
}

#endif