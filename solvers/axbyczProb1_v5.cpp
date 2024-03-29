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

Input:
    A1, B1, C1, A2, B2, C2: Matrices - dim 4x4
    opt: bool
    nstd1, nst2: standard deviation
Output:
    X_final, Y_final, Z_final: Matrices - dim 4x4

In the case of two robotic arms:
     A - robot 1's base to end effector transformation (forward kinematics)
     B - camera to calibration target transformation
     C - robot 2's base to end effector transformation (forward kinematics)
     X - end effector of robot 1 to camera transformation
     Y - robot 1's base to robot 2's base transformation
     Z - end effector of robot 2 to calibration target transformation
*/

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "batchSolveXY.h"
#include "rotError.h"
#include "tranError.h"

void axbyczProb1(const std::vector<Eigen::Matrix4d> &A1,
                 const std::vector<Eigen::Matrix4d> &B1,
                 const std::vector<Eigen::Matrix4d> &C1,
                 const std::vector<Eigen::Matrix4d> &A2,
                 const std::vector<Eigen::Matrix4d> &B2,
                 const std::vector<Eigen::Matrix4d> &C2,
                 bool opt,
                 double nstd1,
                 double nstd2,
                 Eigen::Matrix4d &X_final,
                 Eigen::Matrix4d &Y_final,
                 Eigen::Matrix4d &Z_final) {

    ////   A1 is constant with B1 and C1 free
    Eigen::Matrix4d A1_fixed = A1[0];

    ////   C2 is constant with A2 and B2 free
    Eigen::Matrix4d C2_fixed = C2[0];

    int len = A1.size();

    std::vector<Eigen::Matrix4d> X_g(len), Z_g(len), X(len), Y_temp(len), Z(len);

    Eigen::Matrix4d MeanA, MeanB, MeanC, MeanB1, MeanC1, MeanA2, MeanB2;
    Eigen::Matrix<double, 6, 6> SigA, SigB, SigC;

    std::cout << "A1.size(): " << len << std::endl;

    //// ------ using probability methods ------
    // calculate Z_g : all guesses of Z
    //// ------ Solve for Z -------- //
    // A1 fixed, B1 and C1 free

    batchSolveXY(C1, B1, opt, nstd1, nstd2, Z_g, Y_temp,
                 MeanC1, MeanB1, SigC, SigB);

    std::cout << "Z_g[0]: " << Z_g[0] << std::endl;
    std::cout << "Z_g[9]: " << Z_g[9] << std::endl;
    std::cout << "Y_temp[0]: " << Y_temp[0] << std::endl;
    std::cout << "Y_temp[9]: " << Y_temp[9] << std::endl;

    // Keep the candidates of Z that are SE3
    // Normally there will be four Z \in SE3
    int Z_index = 0;
    for (int i = 0; i < Z_g.size(); ++i) {
        if (Z_g[i].determinant() > 0) {
            Z.push_back(Z_g[i]);
            ++Z_index;
        }
    }

    int s_Z = Z.size();
    std::cout << "s_Z: " << std::endl << s_Z << std::endl;
    std::cout << "Z: " << Z[0] << std::endl;

    //// ------ Solve for X -------- //
    // C2 fixed, A2 and B2 free

    // ------ Calculate B2^-1 -------
    int Num = A2.size();
    std::vector<Eigen::Matrix4d> A2_inv(Num), B2_inv(Num);
    for (int i = 0; i < Num; ++i) {
        A2_inv[i] = A2[i].inverse();
        B2_inv[i] = B2[i].inverse();
    }

    // ------ using probability methods ------
    // calculate X_g : all guesses of X
    batchSolveXY(A2, B2_inv, opt, nstd1, nstd2, X_g, Y_temp,
                 MeanA2, MeanB, SigC, SigB);

    std::cout << "X_g[0]: " << X_g[0] << std::endl;
    std::cout << "X_g[9]: " << X_g[9] << std::endl;
    std::cout << "Y_temp[0]: " << Y_temp[0] << std::endl;
    std::cout << "Y_temp[9]: " << Y_temp[9] << std::endl;

    // Calculate MeanB for computing Y later
    // Note: can be further simplified by using only the distribution function
    batchSolveXY(A2_inv, B2, opt, nstd1, nstd2, X_g, Y_temp,
                 MeanA, MeanB2, SigC, SigB);

    std::cout << "X_g[0]: " << X_g[0] << std::endl;
    std::cout << "X_g[9]: " << X_g[9] << std::endl;
    std::cout << "Y_temp[0]: " << Y_temp[0] << std::endl;
    std::cout << "Y_temp[9]: " << Y_temp[9] << std::endl;

    // Keep the candidates of X that are SE3å
    // Normally there will be four X \in SE3
    int X_index = 0;
    for (int i = 0; i < X_g.size(); ++i) {
        if (X_g[i].determinant() > 0) {
            X.push_back(X_g[i]);
            ++X_index;
        }
    }

    int s_X = X.size();
    std::cout << "s_X: " << std::endl << s_X << std::endl;
    std::cout << "X: " << X[0] << std::endl;

    //// ------ Solve for Y -------- //
    // Compute Y using the mean equations
    std::vector<Eigen::Matrix4d> Y(2*s_X*s_Z);
    for (int i = 0; i < s_X; ++i) {
        for (int j = 0; j < s_Z; ++j) {
            Y[i * s_Z + j] = (A1_fixed * X[i] * MeanB1 * Z[j].inverse()) * MeanC1.inverse();
            Y[i * s_Z + j + s_X * s_Z] = (MeanA2 * X[i] * MeanB2) * Z[j].inverse() * C2_fixed.inverse();
        }
    }

    int s_Y = Y.size();
    std::cout << "s_Y: " << std::endl << s_Y << std::endl;
    std::cout << "Y: " << Y[0] << std::endl;

    //// Find out the optimal (X, Y, Z) that minimizes cost

    Eigen::MatrixXd cost = Eigen::MatrixXd::Zero(s_X, s_Y * s_Z);
    double weight = 1.5; // weight on the translational error of the cost function
    double min_cost = std::numeric_limits<double>::max();
    int min_i = 0, min_j = 0, min_m = 0;

    for (int i = 0; i < s_X; ++i) {
        for (int j = 0; j < s_Z; ++j) {
            for (int m = 0; m < s_Y; ++m) {
                Eigen::MatrixXd left1 = A1_fixed * X[i] * MeanB1;
                Eigen::MatrixXd right1 = Y[m] * MeanC1 * Z[j];

                double diff1 =
                        rotError(left1, right1) + weight * tranError(left1, right1);

                Eigen::MatrixXd left2 = MeanA2 * X[i] * MeanB2;
                Eigen::MatrixXd right2 = Y[m] * C2_fixed * Z[j];
                double diff2 =
                        rotError(left2, right2) + weight * tranError(left2, right2);

                double current_cost = diff1 + diff2;
                cost(i, j * s_Y + m) = current_cost;

                /*std::cout << "X[" << i << "]:\n" << X[i] << std::endl;
                std::cout << "Y[" << m << "]:\n" << Y[m] << std::endl;
                std::cout << "Z[" << j << "]:\n" << Z[j] << std::endl;*/

                if (current_cost < min_cost) {
                    min_cost = current_cost;
                    min_i = i;
                    min_j = j;
                    min_m = m;
                }
            }
        }
    }

    std::cout << "Minimum cost: " << min_cost << std::endl;
    std::cout << "Indices (min_i, min_j, min_m): (" << min_i << ", " << min_j << ", " << min_m << ")" << std::endl;

    std::cout << "X_final: " << std::endl << X[1] << std::endl;
    std::cout << "Y_final: " << std::endl << Y[1] << std::endl;
    std::cout << "Z_final: " << std::endl << Z[1] << std::endl;
    std::cout << "X_final: " << std::endl << X[2] << std::endl;
    std::cout << "Y_final: " << std::endl << Y[2] << std::endl;
    std::cout << "Z_final: " << std::endl << Z[2] << std::endl;

    //// recover the X,Y,Z that minimizes cost
    X_final = X[min_i];
    Z_final = Z[min_j];
    Y_final = Y[min_m];
}

int main() {
    int num = 10;
    srand(12345);

    std::vector<Eigen::Matrix4d> A1(num), B1(num), C1(num),
                                    A2(num), B2(num), C2(num);

    for (int i = 0; i < num; ++i){
        A1[i] = Eigen::Matrix4d::Random();
        B1[i] = Eigen::Matrix4d::Random();
        C1[i] = Eigen::Matrix4d::Random();
        A2[i] = Eigen::Matrix4d::Random();
        B2[i] = Eigen::Matrix4d::Random();
        C2[i] = Eigen::Matrix4d::Random();

        A1[i].row(3) << 0, 0, 0, 1;
        B1[i].row(3) << 0, 0, 0, 1;
        C1[i].row(3) << 0, 0, 0, 1;
        A2[i].row(3) << 0, 0, 0, 1;
        B2[i].row(3) << 0, 0, 0, 1;
        C2[i].row(3) << 0, 0, 0, 1;
    }

    bool opt = true;
    double nstd1 = 0.01;
    double nstd2 = 0.01;

    Eigen::Matrix4d X_final;
    Eigen::Matrix4d Y_final;
    Eigen::Matrix4d Z_final;

    axbyczProb1(A1,B1,C1,A2,B2,C2,opt,nstd1,nstd2,X_final,Y_final,Z_final);

    std::cout << "Build successful? - YES" <<std::endl;

    std::cout << "X_final: " << std::endl << X_final << std::endl;
    std::cout << "Y_final: " << std::endl << Y_final << std::endl;
    std::cout << "Z_final: " << std::endl << Z_final << std::endl;

    return 0;
}

/*
 * Problem - doesn't seem to work when:
 *  opt = false
 *  nstd1 = 0
 *  ntdd2 = 0
 *  It's working for the above conditions now
 *  New output:
 X_final_:
 -0.373014     -0.789  -0.488201     9.5636
-0.0133807   0.530697  -0.847456   -29.8597
  0.927729  -0.309581  -0.208515   -10.6073
         0          0          0          1

Y_final_:
nan nan nan nan
nan nan nan nan
nan nan nan nan
nan nan nan nan

Z_final_:
-0.833256  0.503671 -0.228036  -16.2424
 0.552888  0.759081 -0.343673  -333.133
        0 -0.412446 -0.910982   281.394
        0         0         0         1


Output:
len: 8
works till here - solve for Z? - YES
works till here - solve for Z and X? - YES
Y:
-2.58501 0.864086  1.20978 -1.80958
 2.31146 0.523348 0.238907 0.646305
-1.47282 -0.21759  2.77298 -1.66964
-5.26496 -2.34186 0.531769 -1.83458
works till here - solve for Z, X and Y? - YES
works till here - optimal cost? - YES
X_final_:
          1           0           0           0
          0   -0.894427   -0.447214 4.80553e-49
          0    0.447214   -0.894427 9.61107e-49
          0           0           0           1
Y_final_:
 -0.808833  -0.843726   -2.20955    1.53485
  0.458386 -0.0436236   0.150767  -0.366332
  -1.15638   -2.61803   -1.09865    1.13648
 -0.179155  -0.866537   0.647998    1.02978
Z_final_:
         -1           0           0 3.38334e-32
          0    0.786109   -0.618088 4.90051e-32
          0   -0.618088   -0.786109  1.0926e-31
          0           0           0           1
works till here - recover X,Y,Z final? - YES
Build successful? - YES

 std::cout << "s_X size: " << s_X << std::endl;
    std::cout << "s_Y size: " << s_Y << std::endl;
    std::cout << "s_Z size: " << s_Z << std::endl;
    std::cout << "cost size: " << cost.size() << std::endl;
    std::cout << "cost rows: " << cost.rows() << std::endl;
    std::cout << "cost columns: " << cost.cols() << std::endl;
    std::cout << "cos(0,0) = " << cost(0, 0) << std::endl;

 */