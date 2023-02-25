/*
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
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <random>
#include <expm.h>
#include <se3Vec.h>
#include <so3Vec.h>

using namespace Eigen;

std::vector<MatrixXd> sensorNoise(const std::vector<MatrixXd> &g, const MatrixXd &gmean, const double &std, const int &model)
{
    std::vector<MatrixXd> g_noise(g.size());

    switch (model)
    {
        case 1:
        {
            // Independently from Normal Distribution
            for (int i = 0; i < g.size(); i++)
            {
                VectorXd temp = VectorXd::Random(3);

                VectorXd noise_old1 = gmean + VectorXd::Zero(6);
                noise_old1.tail(3) = std * VectorXd::Random(3);

                VectorXd noise_old2 = gmean + VectorXd::Zero(6);
                noise_old2.head(3) = std * (temp / temp.norm());
                
                MatrixXd g_temp = g[i] * expm(se3Vec(noise_old1)) * expm(se3Vec(noise_old2));
                g_noise[i] = g_temp;
            }
            break;
        }

        case 2:
        {
            // Coupling Matrix
            MatrixXd C = MatrixXd::Zero(g.size() * 6, g.size() * 6);
            for (int i = 0; i < g.size(); i++)
            {
                if ((i - 3) >= 0)
                    C.block(i * 6, (i - 3) * 6, 6, 6) = MatrixXd::Identity(6, 6) * 0.25;
                if ((i - 2) >= 0)
                    C.block(i * 6, (i - 2) * 6, 6, 6) = MatrixXd::Identity(6, 6) * 0.5;
                C.block(i * 6, i * 6, 6, 6) = MatrixXd::Identity(6, 6);
                if ((i + 1) <= g.size())
                    C.block(i * 6, (i + 1) * 6, 6, 6) = MatrixXd::Identity(6, 6) * 0.5;
                if ((i + 2) <= g.size())
                    C.block(i * 6, (i + 2) * 6, 6, 6) = MatrixXd::Identity(6, 6) * 0.25;
            }

            VectorXd noise_old = std * VectorXd::Random(g.size() * 6) + gmean;
            VectorXd noise_new = C * noise_old;

            for (int i = 0; i < g.size(); i++)
            {
                MatrixXd g_temp = g[i] * expm(se3Vec(noise_new.segment(i * 6, 6)));
                g_noise[i] = g_temp;
            }
            break;
        }

        case 3:
        {
            // Wiener Process
            MatrixXd noise_old = std * MatrixXd::Random(6, g.size()) + gmean;
            MatrixXd noise_new = MatrixXd::Zero(6, g.size());
            
            for (int i = 0; i < g.size(); i++)
            {
                noise_new.col(i).head(6) = noise_old.leftCols(i + 1).rowwise().sum() / sqrt(i + 1);
                MatrixXd g_temp = g[i] * expm(se3Vec(noise_new.col(i).head(6)));
            }
        }

        case 4:
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
        }    

        case 5:
        {
            //-----------outlier noise-----------------------------------
            for (int i = 0; i < g_noise.size(); i++) {
                Matrix4d noise_old1, noise_old2;

                if ((i + 1) % 10 == 1) {
                    Vector3d temp = Vector3d::Random();
                    noise_old1 << Matrix3d::Zero(), 0.1*Matrix3d::Zero(), temp + gmean, Vector4d::UnitW();
                    noise_old2 << 0.01*(temp / temp.norm()), Matrix3d::Zero(), gmean, Vector4d::UnitW();
                    g_noise[i] = g[i]*expm(se3Vec(noise_old1))*expm(se3Vec(noise_old2));
                }
                else {
                    noise_old1 << std*Matrix3d::Identity(), Matrix3d::Zero(), gmean, Vector4d::UnitW();
                    noise_old2 << Matrix3d::Zero(), std*Matrix3d::Identity(), gmean, Vector4d::UnitW();
                    g_noise[i] = g[i]*expm(se3Vec(noise_old1))*expm(se3Vec(noise_old2));
                }
            }
        }

        default:
            std::cout<<"Error"<<std::endl;
    }
}