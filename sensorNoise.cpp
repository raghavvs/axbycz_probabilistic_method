#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <random>
#include "se3Vec.h" // assuming that the se3_vec() function has been implemented in a separate file

using namespace std;

// Function to add Gaussian noise to a g in SE(3)
// model = 1: Add noise independently from Normal Distribution
// g: input SE(3) matrices to add noise to
// gmean: mean value of the noise
// std: standard deviation of the noise
// g_noise: output SE(3) matrices with added noise
void sensorNoise(MatrixXd& g, const VectorXd& gmean, double std, int model, MatrixXd& g_noise)
{
    random_device rd;  // obtain a random seed from the hardware
    mt19937 gen(rd()); // seed the generator
    normal_distribution<> d(0, std); // define the normal distribution with mean 0 and standard deviation std
    
    switch (model) {
        case 1:
            // Add noise independently from Normal Distribution
            for (int i = 0; i < g.cols(); i++) {
                VectorXd temp = VectorXd::Random(3); // generate a random vector of size 3
                VectorXd noise_old1 = VectorXd::Zero(6); // initialize the noise_old1 vector to [0, 0, 0, 0, 0, 0]
                VectorXd noise_old2 = VectorXd::Zero(6); // initialize the noise_old2 vector to [0, 0, 0, 0, 0, 0]
                for (int j = 0; j < 3; j++) {
                    noise_old1(j+3) = d(gen) + gmean(j+3); // add noise to the translational component of the SE(3) matrix
                    noise_old2(j) = d(gen)*temp(j)/temp.norm() + gmean(j); // add noise to the rotational component of the SE(3) matrix
                }
                g_noise.col(i) = g.col(i) * expm(se3_vec(noise_old1)) * expm(se3_vec(noise_old2)); // add noise to the SE(3) matrix
            }
            break;
            
        default:
            cerr << "Invalid model type. Please use model = 1." << endl;
            break;
    }
}
