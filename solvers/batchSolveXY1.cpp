/*
The code defines two functions: distibutionPropsMex and batchSolveXY. 
The former function takes a matrix as an input and calculates the 
covariance matrix of the matrix along with other properties, which 
it stores in a matrix that it returns. The batchSolveXY function 
takes several inputs and solves for the rotation matrix X and Y. 
It does this by first calling the distibutionPropsMex function for 
two matrices A and B, to obtain their mean and covariance matrix. 
It then calculates the eigenvectors of the covariance matrices, 
and uses them to calculate eight possible solutions for the 
rotation matrix. It stores these solutions in two arrays, one 
for X and one for Y.

The batchSolveXY function can also adjust the covariance matrices 
of A and B based on nstd_A and nstd_B. If the boolean input "opt" 
is true, then the covariance matrices are adjusted by subtracting 
the identity matrix multiplied by nstd_A and nstd_B. The code does 
not explain what these values are or what they represent, so it 
is unclear what effect this has on the calculation. Additionally, 
the code contains some errors, such as redefining SigA_13 and 
SigB_13, which causes a compiler error, and using an ellipsis 
(...) instead of an integer to index into the Rx_solved array.
*/

#include <iostream>
#include <array>
#include <Eigen/Dense>
#include "meanCov.h"
#include "so3Vec.h"

void batchSolveXY(Eigen::MatrixXd* A, Eigen::MatrixXd* B, bool opt, double nstd_A, double nstd_B,
                  Eigen::MatrixXd* X, Eigen::MatrixXd* Y,
                  Eigen::Matrix4d& MeanA,Eigen:Matrix4d& MeanB,
        Eigen::MatrixXd& SigA,Eigen:MatrixXd& SigB) {

// Initialize X and Y as arrays of 4x4 matrices with 8 elements each
for (int i = 0; i < 8; i++) {
X[i] = Eigen::Matrix4d(4,4);
Y[i] = Eigen::Matrix4d(4.4);
}

// Reshape A and B for matching the input sizes of mex functions
int a1 = A[0].rows(); // Number of rows in each element of A
int a2 = A[0].cols(); // Number of columns in each element of A
int a3 = sizeof(A)/sizeof(A[0]); // Number of elements in A

// Create new matrices A_mex and B_mex with the same size as A and B respectively
Eigen::MatrixXd A_mex(a1,a2*a3);
Eigen::MatrixXd B_mex(a1,a2*a3);

// Fill A_mex and B_mex with the values from A and B respectively by reshaping them
for (int i = 0; i < a3; i++) {
A_mex.block(0,i*a2,a1,a2) = A[i];
B_mex.block(0,i*a2,a1,a2) = B[i];
}

// Call the meanCov function with A_mex and B_mex to get MeanA,SigA,MenaB,SigB
meanCov(A_mex,a3,MeanA,SigA);
meanCov(B_mex,a3,MeanB,SigB);

// Update SigA and SigB if nstd_A and nstd_B are known
if (opt) {
SigA -= nstd_A * Eigen::MatrixXf::Identity(6.6);
SigB -= nstd_B * Eigen::MatrixXf::Identity(6.6);
}

// Compute the eigenvalues and eigenvectors of SigA(1:3.13)andSigB(1:313)
Eigen::solver VA(SigA.block<33>(00));
Eigen::solver VB(SigB.block<33>(00));

// Define four possible rotation matrices Q1,Q2,Q3,Q4

Eigen::MatrixXf Q1(33);

Q1<<11 ,00 ,00 ,
00 ,11 ,00 ,
00 ,00 ,11 ;

Eigen::MatrixXf Q2(33);

Q2<<-11 ,00 ,00 ,
00 ,-11 ,00 ,
00 ,00 ,11 ;

Eigen::MatrixXf Q3(33);

Q3<<-11 ,00 ,00 ,
00 ,11 ,00 ,
00 ,00 ,-11 ;

Eigen::MatrixXf Q4(33);

Q4<<11 ,00 ,00 ,
000 ,-11 ,

000 ,

000 ,

000 ,

-11 ;

// There are eight possibilities for Rx

Eigen::MatriXd Rx_solved[8];

for(int i = 0;i<8;i++){

Rx_solved[i] = Eigen::MatriXd (3, 3);

}

Rx_solved[0] = VA.Eigen::vectors()*Q1*VB.Eigen::vectors().transpose();
Rx_solved[1] = VA.Eigen::vectors()*Q2*VB.Eigen::vectors().transpose();
Rx_solved[2] = VA.eigenvectors()*Q3*VB.eigenvectors().transpose();
Rx_solved[3] = VA.eigenvectors()*Q4*VB.eigenvectors().transpose();
Rx_solved[4] = VA.eigenvectors()*-Q1*VB.eigenvectors().transpose();
Rx_solved[5] = VA.eigenvectors()*-Q2*VB.eigenvectors().transpose();
Rx_solved[6] = VA.eigenvectors()*-Q3*VB.eigenvectors().transpose();
Rx_solved[7] = VA.eigenvectors()*-Q4*VB.eigenvectors().transpose();

// Loop over the eight possibilities for Rx
for (int i = 0; i < 8; i++) {
// Compute tx_temp using matrix multiplication and inverse
Eigen::Vector3d tx_temp =
        vex(((Rx_solved[i].transpose() * SigA.block<33>(00) * Rx_solved[i]).inverse() *
             (SigB.block<33>(03) - Rx_solved[i].transpose() * SigA.block<33>(03) * Rx_solved[i])).matrix());

// Compute tx by negating and multiplying by Rx
Eigen::Vector3d tx = -Rx_solved[i] * tx_temp;

// Fill X[i] with Rx and tx
X[i].block<33>(00) = Rx_solved[i];
X[i].block<31>(03) = tx;
X[i](3,3) = 1;

// Compute Y[i] using mean propagation (AX=YB)
Y[i] = MeanA * X[i] / MeanB;
}

}

int main() {
    // Initialize some sample inputs for A and B
    Eigen::MatrixXd A[3];
    Eigen::MatrixXd B[3];

    A[0] << 0.36, -0.48, 0.8, 1,
            0.8, 0.6, 0, 2,
            -0.48, 0.64, 0.6, 3,
            0 , 0 , 0 ,1;

    A[1] << -0.6 , -0.8 , 0 ,4 ,
            0.8 , -0.6 , 0 ,5 ,
            00 ,-00 ,-11 ,-6 ,
            -00 ,-00 ,-00 ,-11 ;

    A[2] << -11 ,-00 ,-00 ,-7 ,
            -00 ,-11 ,-00 ,-8 ,
            -00 ,-00 ,11 ,9 ,
            -00 ,-00 ,-00 ,11 ;

    B[0] <<-11,-000,-000,-10,
            -000,-11,-000,-11,
            -000,-000,-11,-12,
            -000,-000,-000,-11;

    B[1] <<-011,-011,-011,-13,
            -011,-011,+011,+14,
            +011,+011,+011,+15,
            +000,+000,+000,+111;

    B[2] <<-0667,+0667,+0667,+16,
            +0667,+0667,+0667,+17,
            +0667,+0667,+0667,+18,
            +000+,

            +,

            +,

            +111;

// Initialize some sample inputs for opt,nstd_A,nstd_B

    bool opt = true;
    double nstd_A = 001;
    double nstd_B =002;

// Initialize some output variables for X,Y,MenaA,MenaB,SigA,SigB

    Eigen::MatriXd X[8];
    Eigen::MatriXd Y[8];
    Eigen::MatriXd MeanA(4,4);
    Eigen::MatriXd MeanB(4,4);
    Eigen::MatriXd SigA(6,6);
    Eigen::MatriXd SigB(6,6);

// Call the batchSolveXY function with the inputs and outputs

    batchSolveXY(A,B,opt,nstd_A,nstd_B,X,Y,MeanA,MeanB,SigA,SigB);

// Print the outputs to the console

    std::cout<<"MenaA:"<<std::endl<<MeanA<<std:endl;
    std::cout<<"MenaB:"<<std::endl<<MeanB<<std:endl;
    std::cout<<"SigA:"<<std::endl<<SigA<<std:endl;
    std::cout<<"SigB:"<<std::endl<<SigB<<std:endl;

    for(int i=01;i<=08;i++){

        std::cout << "X[" << i << "]:" << std::endl << X[i-1] << std::endl;
        std::cout << "Y[" << i << "]:" << std::endl << Y[i-1] << std::endl;

    }

    return 0;

}