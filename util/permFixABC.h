/*
Generate partially permuted data triples [A1_perm, B1_perm, C1_perm]
 Input:
       M in 4 x 4 x 1
       N in 4 x 4 x n
       P in 4 x 4 x n
       r : scrambling rate
 Output:
       M_perm in 4 x 4 x n (n copies of M)
       N_perm in 4 x 4 x n (permuated N)
       P_perm in 4 x 4 x n (samme as P)
*/

/*
DESCRIPTION:

This program defines a function called permFixABC that takes as input three matrices 
M, N, and P, and a scalar r. It generates partially permuted data triples [A1_perm, B1_perm, C1_perm]. 
M is a 4x4x1 matrix, N is a 4x4xn matrix, P is a 4x4xn matrix, and r is a scrambling rate. 
The output consists of three matrices, M_perm which is a copy of M, N_perm which is a 
permuted version of N, and P_perm which is the same as P. The permuted version of N is 
obtained by calling the function scrambleData with N and r. The function permFixABC uses 
Eigen matrices and blocks to resize and fill the matrices M_perm.
*/

#ifndef PERMFIXABC_H
#define PERMFIXABC_H


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <scrambleData.h>

void permFixABC(const Eigen::MatrixXd& M, const Eigen::MatrixXd& N, const Eigen::MatrixXd& P, 
                double r, Eigen::MatrixXd& M_perm, Eigen::MatrixXd& N_perm, Eigen::MatrixXd& P_perm) {

    int n = N.cols();
    M_perm.resize(M.rows(), M.cols() * n);
    for (int i = 0; i < n; i++) {
        M_perm.block(0, i * M.cols(), M.rows(), M.cols()) = M;
    }
    N_perm = scrambleData(N, r);
    P_perm = P;
}

#endif