#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>

using namespace Eigen;
using namespace std;

double metric(vector<Matrix4d>& A, vector<Matrix4d>& B, vector<Matrix4d>& C, 
              Matrix4d& X, Matrix4d& Y, Matrix4d& Z) {
    double diff = 0.0;
    int N = 0;
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].cols(); ++j) {
            diff += (A[i].block<3, 3>(0, 0) * X * B[i].block<3, 3>(0, 0) - Y * C[i].block<3, 3>(0, 0) * Z).norm();
            N++;
        }
    }
    return diff / N;
}
