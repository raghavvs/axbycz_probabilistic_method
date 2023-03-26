/*
DESCRIPTION:

The program defines a function scrambleData that takes a 3D Eigen::MatrixXd M 
as input along with a scrambling rate s_rate. It then generates a random 
permutation of the third dimension of M with a probability proportional to 
s_rate and returns a partially permuted copy of M as output. The program 
uses the C++ standard library chrono and random to generate the random permutation.
*/

#ifndef SCRAMBLEDATA_H
#define SCRAMBLEDATA_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <random>
#include <algorithm>
#include <chrono>

std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
std::uniform_real_distribution<double> distribution(0.0, 1.0);

Eigen::MatrixXd scrambleData(Eigen::MatrixXd M, double s_rate) {
    int n = M.cols();
    std::vector<int> M_index(n);
    for (int i = 0; i < n; i++) {
        M_index[i] = i;
    }
    for (int i = 0; i < n; i++) {
        if (distribution(generator) <= 0.01 * s_rate) {
            int index = rand() % n;
            std::swap(M_index[i], M_index[index]);
        }
    }
    Eigen::MatrixXd M_perm = Eigen::MatrixXd::Zero(M.rows(), M.cols());
    for (int i = 0; i < n; i++) {
        M_perm.col(i) = M.col(M_index[i]);
    }
    return M_perm;
}

#endif