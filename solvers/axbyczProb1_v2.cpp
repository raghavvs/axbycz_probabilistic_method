/*
DESCRIPTION:

The provided code implements a set of functions that solve a robotics problem 
involving the transformation matrices of multiple coordinate frames. Specifically, 
the functions solve for the transformations between three coordinate frames 
(A, B, and C) given the transformations between A and B and between A and C. 
The functions use the Eigen library to perform matrix operations such as inversion 
and SVD decomposition. The main function (axbyczProb1) calls the other two functions 
(batchSolveXY and randSE3) to generate a set of random transformations and iteratively 
select those that satisfy certain constraints, in order to estimate the desired transformations.
*/

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <batchSolveXY.h>

using namespace Eigen;

void axbyczProb1(Matrix4d A1, Matrix4d B1, Matrix4d C1, Matrix4d A2, Matrix4d B2, Matrix4d C2,
                  int opt, double nstd1, double nstd2, Matrix4d& X_final, Matrix4d& Y_final, Matrix4d& Z_final)
{
    Matrix4d A1_1 = A1, C2_1 = C2;

    Matrix4d Z_g[4];
    Matrix4d MeanC1, MeanB1;
    int s_Z = 0;

    // calculate Z_g : all guesses of Z
    // Normally, there will be four Z ∈ SE(3)
    // Keep the candidates of Z that are SE(3)
    for (int i = 0; i < 4; i++)
    {
        // Note that the MATLAB function `batchSolveXY` is not defined here.
        // You'll need to provide an implementation of this function.
        // This function should take two 4x4 matrices and return a 4x4 matrix.
        Z_g[i] = batchSolveXY(C1, B1, opt, nstd1, nstd2);
        if (Z_g[i].determinant() > 0)
        {
            s_Z++;
        }
    }

    Matrix4d X_g[4];
    Matrix4d MeanA2, MeanB2;
    int s_X = 0;

    // Calculate B2^-1
    int Num = A2.cols();
    Matrix4d A2_inv[Num], B2_inv[Num];
    for (int i = 0; i < Num; i++)
    {
        A2_inv[i] = A2[i].inverse();
        B2_inv[i] = B2[i].inverse();
    }

    // calculate X_g : all guesses of X
    // Keep the candidates of X that are SE(3)
    for (int i = 0; i < 4; i++)
    {
        X_g[i] = batchSolveXY(A2, B2_inv, opt, nstd1, nstd2);
        if (X_g[i].determinant() > 0)
        {
            s_X++;
        }
    }

    // Compute Y using the mean equations
    Matrix4d Y[2 * s_X * s_Z];
    int k = 0;
    for (int i = 0; i < s_X; i++)
    {
        for (int j = 0; j < s_Z; j++)
        {
            // There are at least four mean equations to choose from to compute
            // Y. It will be interesting to see how each choice of the mean
            // equations can affect the result
            Y[k] = (A1_1.array() * X_g[i].array() * MeanB1.array() / Z_g[j].array()).matrix() / MeanC1.array();
            Y[k + s_X * s_Z] = (MeanA2.array() * X_g[i].array() * MeanB2.array() / Z_g[j].array()).matrix() / C2_1.array();
            k++;
        }
    }

    // Find out the optimal (X, Y, Z) that minimizes cost
    int s_Y = 2 * s_X * s_Z;
    MatrixXd cost(s_X, s_Y * s_Z);
    double weight = 1.0 / (2.0 * sigma * sigma);

    for (int x = 0; x < s_X; x++) {
        for (int yz = 0; yz < s_Y * s_Z; yz++) {
            int y = yz % s_Y;
            int z = yz / s_Y;
            double diff = (data(x, y, z) - mu(x, y, z));
            cost(x, yz) = weight * diff * diff;
        }
    }

    MatrixXd D(s_X, s_Y * s_Z);
    D.setConstant(0.0);

    for (int x = 1; x < s_X; x++) {
        for (int yz = 0; yz < s_Y * s_Z; yz++) {
            int y = yz % s_Y;
            int z = yz / s_Y;

            double min_val = D(x - 1, yz);
            if (y > 0 && z > 0) {
                min_val = std::min(min_val, D(x - 1, (y - 1) * s_Z + (z - 1)));
            }
            if (y > 0) {
                min_val = std::min(min_val, D(x - 1, (y - 1) * s_Z + z));
            }
            if (z > 0) {
                min_val = std::min(min_val, D(x - 1, y * s_Z + (z - 1)));
            }

            D(x, yz) = cost(x, yz) + min_val;
        }
    }

    // Backtracking to find the optimal path
    int best_yz = 0;
    for (int yz = 1; yz < s_Y * s_Z; yz++) {
        if (D(s_X - 1, yz) < D(s_X - 1, best_yz)) {
            best_yz = yz;
        }
    }

    VectorXi seam_indices(s_X);
    seam_indices(s_X - 1) = best_yz % s_Y;
    for (int x = s_X - 2; x >= 0; x--) {
        int y = seam_indices(x + 1);
        int z = best_yz / s_Y;
        double min_val = D(x, y * s_Z + z);
        seam_indices(x) = y;
        if (y > 0 && z > 0) {
            if (D(x, (y - 1) * s_Z + (z - 1)) < min_val) {
                min_val = D(x, (y - 1) * s_Z + (z - 1));
                seam_indices(x) = y - 1;
            }
        }
        if (y > 0) {
            if (D(x, (y - 1) * s_Z + z) < min_val) {
                min_val = D(x, (y - 1) * s_Z + z);
                seam_indices(x) = y - 1;
            }
        }
        if (z > 0) {
            if (D(x, y * s_Z + (z - 1)) < min_val) {
                min_val = D(x, y * s_Z + (z - 1));
                seam_indices(x) = y;
            }
        }
        best_yz = seam_indices(x) * s_Y + z;
    }

    return seam_indices;

}

int main() {
    // Load the input data
    MatrixXd
    input_data;
    ifstream input_file("input.txt");
    input_file >> input_data;
    // Get the dimensions of the input data
    int n = input_data.rows();
    int m = input_data.cols();

    // Compute the optimal (X, Y, Z) that minimizes cost
    int s_X = round(pow(n, 1.0 / 3.0));
    int s_Z = round(pow(m, 1.0 / 3.0));
    int s_Y = 2 * s_X * s_Z;
    MatrixXd cost(s_X, s_Y * s_Z);
    double weight = 1.0 / (s_X * s_Y * s_Z);

    for (int i = 0; i < s_X; i++) {
        for (int j = 0; j < s_Y * s_Z; j++) {
            cost(i, j) = weight * compute_cost(input_data, s_X, s_Y, s_Z, i, j);
        }
    }

    VectorXd X, Y, Z;
    tie(X, Y, Z) = optimize_cost(cost);

    // Compute the tensor approximation
    MatrixXd tensor_approximation = compute_tensor_approximation(input_data, X, Y, Z);

    // Save the tensor approximation to a file
    ofstream output_file("output.txt");
    output_file << tensor_approximation;

    return 0;
}
