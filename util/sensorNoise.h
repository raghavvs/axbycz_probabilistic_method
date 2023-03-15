/*
DESCRIPTION:

Function to add Gaussian noise to a g in SE(3)
model = 1: Add noise independently from Normal Distribution
g: input SE(3) matrices to add noise to
gmean: mean value of the noise
std: standard deviation of the noise
g_noise: output SE(3) matrices with added noise

This code includes a C++ function called sensorNoise which takes four input arguments: 
a vector of matrices g, a matrix gmean, a double std, and an integer model. 
The output is a vector of matrices g_noise.

The function applies different noise models to the input matrices g based on 
the value of the input integer model. If model is 1, then the function applies 
independent Gaussian noise to each matrix. If model is 2, then the function 
applies a coupling matrix to introduce correlation between the noise applied 
to each matrix. If model is 3, then the function applies Wiener process noise 
to each matrix. If model is 4, then the function applies Plücker nudge noise 
to each matrix.

The function returns the input matrices with the applied noise as a vector of matrices.
New argument - length of the matrix array

ONLY CASE 1 - WORKS - at the moment
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include <random>
#include "se3Vec.h"
#include "so3Vec.h"

Eigen::Matrix4d* sensorNoise(const Eigen::Matrix4d* g, int len, const Eigen::MatrixXd gmean, double sd, int model)
    {
        // Declare g_noise as an array of matrices and allocate memory for it
        Eigen::Matrix4d *g_noise = new Eigen::Matrix4d[len];

        switch (model) {
            case 1: {
                Eigen::Vector3d temp = Eigen::Vector3d::Random();
                Eigen::VectorXd noise_old1(6);
                Eigen::VectorXd noise_old2(6);

                // Independently from Normal Distribution
                noise_old1.segment(0, 3) = Eigen::Vector3d::Zero();
                noise_old1.segment(3, 3) = sd * Eigen::Vector3d::Random();
                noise_old1 += gmean;

                noise_old2.segment(0, 3) = sd * temp.normalized() * temp.norm();
                noise_old2.segment(3, 3) = Eigen::Vector3d::Zero();
                noise_old2 += gmean;

                Eigen::Matrix4d exp1 = (se3Vec(noise_old1)).exp();
                Eigen::Matrix4d exp2 = (se3Vec(noise_old2)).exp();
                Eigen::Matrix4d prod = exp1 * exp2;

                for (int i = 0; i < len; i++) {
                    g_noise[i] = g[i] * prod;
                }

                break;
            }

                /*case 2:
                {
                    // Coupling Matrix
                    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(g.size() * 6, g.size() * 6);
                    for (int i = 0; i < g.size(); i++)
                    {
                        if ((i - 3) >= 0)
                            C.block(i * 6, (i - 3) * 6, 6, 6) = Eigen::MatrixXd::Identity(6, 6) * 0.25;
                        if ((i - 2) >= 0)
                            C.block(i * 6, (i - 2) * 6, 6, 6) = Eigen::MatrixXd::Identity(6, 6) * 0.5;
                        C.block(i * 6, i * 6, 6, 6) = Eigen::MatrixXd::Identity(6, 6);
                        if ((i + 1) < g.size())
                            C.block(i * 6, (i + 1) * 6, 6, 6) = Eigen::MatrixXd::Identity(6, 6) * 0.5;
                        if ((i + 2) < g.size())
                            C.block(i * 6, (i + 2) * 6, 6, 6) = Eigen::MatrixXd::Identity(6, 6) * 0.25;
                    }

                    Eigen::VectorXd noise_old = sd * Eigen::VectorXd::Random(g.size() * 6) + gmean;
                    Eigen::VectorXd noise_new = C * noise_old;

                    for (int i = 0; i < g.size(); i++)
                    {
                        Eigen::MatrixXd g_temp = g[i] * (se3Vec(noise_new.segment(i * 6, 6))).exp();
                        g_noise[i] = g_temp;
                    }
                    break;
                }


                case 3:
                {
                    // Wiener Process
                    Eigen::MatrixXd noise_old = sd * Eigen::MatrixXd::Random(g.size(), 6) + gmean;
                    Eigen::MatrixXd noise_new = Eigen::MatrixXd::Zero(6, g.size());

                    for (int i = 0; i < g.size(); i++)
                    {
                        noise_new.col(i).head(6) = noise_old.leftCols(i + 1).rowwise().sum() / sqrt(i + 1);
                        Eigen::MatrixXd g_temp = g[i] * (se3Vec(noise_new.col(i).head(6))).exp();
                        g_noise[i] = g_temp;
                    }
                    break;
                }

                case 4:
                {
                    // Outlier noise
                    for (int i = 0; i < g_noise.size(); i++) {
                        Eigen::Matrix4d noise_old1, noise_old2;

                        if ((i + 1) % 10 == 1) {
                            Eigen::Vector3d temp = Eigen::Vector3d::Random();
                            noise_old1 << Eigen::Matrix3d::Zero(), 0.1*Eigen::Matrix3d::Zero(), temp + gmean, Eigen::Vector4d::UnitW();
                            noise_old2 << 0.01*(temp / temp.norm()), Eigen::Matrix3d::Zero(), gmean, Eigen::Vector4d::UnitW();
                            g_noise[i] = g[i] * (se3Vec(noise_old1)).exp() * (se3Vec(noise_old2)).exp();
                        }
                        else {
                            noise_old1 << sd * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(), gmean, Eigen::Vector4d::UnitW();
                            noise_old2 << Eigen::Matrix3d::Zero(), sd * Eigen::Matrix3d::Identity(), gmean, Eigen::Vector4d::UnitW();
                            g_noise[i] = g[i] * (se3Vec(noise_old1)).exp() * (se3Vec(noise_old2)).exp();
                        }
                    }
                    break;
                }*/

            default:
                std::cout << "Error" << std::endl;
        }
        return g_noise;
    }

// TEST CASE

/* int main()
{
    // Create example input data
    Eigen::MatrixXd g1 = Eigen::MatrixXd::Identity(4, 4);
    Eigen::MatrixXd g2 = Eigen::MatrixXd::Identity(4, 4);
    std::vector<Eigen::MatrixXd> g = {g1, g2};
    Eigen::MatrixXd gmean = Eigen::MatrixXd::Zero(6, 1);
    double std = 0.1;
    int model = 1;

    // Call the sensorNoise function
    std::vector<Eigen::MatrixXd> g_noise = sensorNoise(g, gmean, std, model);

    // Print the original and noisy data
    std::cout << "Original data:\n";
    for (const auto& g_matrix : g)
    {
        std::cout << g_matrix << "\n\n";
    }

    std::cout << "Noisy data:\n";
    for (const auto& g_matrix : g_noise)
    {
        std::cout << g_matrix << "\n\n";
    }

    return 0;
} */





/*         case 4: Plucker Nudge - errors - needs debugging - at the moment - ignoring this case
        {   
            //-----------plucker nudge------------------------------------
            for (int i = 0; i < g.size(); i++) {

                Matrix3d Ng, Ng2;
                Vector3d pg, pg2, n, u;
                double thetag, thetag2, dg, dg2; 

                double param_extract(g[i].col(i), thetag, Ng, dg, pg);

                thetag2 = thetag + std*randn();
                Ng2 = so3Vec(so3Vec(Ng) + (Vector3d() << std*randn(), std*randn(), std*randn()).finished().normalized()*std);
                Matrix4d Rg2 = (Matrix4d() << expm(thetag2*Ng2), Vector3d::Zero(), 0, 0, 0, 1).finished();
                Vector3d tg2 = g[0](0, 4 * i) + (Vector3d() << std*randn(), std*randn(), std*randn()).finished();
                dg2 = tg2.dot(so3Vec(Ng2));

                n = so3Vec(Ng2);
                u = (Vector3d() << -n(1), n(0), 0).finished().normalized();
                Vector2d c = (Vector2d() << std::cos(thetag2) - 1, std::sin(thetag2)).finished().inverse() * (Vector2d() << tg2.dot(u), tg2.dot(n.cross(u))).finished();

                pg2 = c(0)*u + c(1)*n.cross(u);

                Matrix4d g_noise_i;
                g_noise_i << expm(thetag2*Ng2), (Matrix3d::Identity() - expm(thetag2*Ng2))*pg2 + dg2*so3Vec(Ng2), 0, 0, 0, 1;
                g_noise[i] = g_noise_i;
            }
        }  */