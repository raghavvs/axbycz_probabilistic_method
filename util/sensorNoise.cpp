/*
Function to add Gaussian noise to a g in SE(3)
model = 1: Add noise independently from Normal Distribution
g: input SE(3) matrices to add noise to
gmean: mean value of the noise
std: standard deviation of the noise
g_noise: output SE(3) matrices with added noise
*/

// INCOMPLETE

#include <cmath>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
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
                
                MatrixXd g_temp = g[i] * expm(se3_vec(noise_old1)) * expm(se3_vec(noise_old2));
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
                MatrixXd g_temp = g[i] * expm(se3_vec(noise_new.segment(i * 6, 6)));
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
                MatrixXd g_temp = g[i] * expm(se3_vec

