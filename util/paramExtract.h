//
// Created by Raghavendra N S on 6/7/23.
//

#ifndef PARAMEXTRACT_H
#define PARAMEXTRACT_H

#include <Eigen/Dense>

void paramExtract(double& theta,
                   Eigen::MatrixXd& N,
                   double& d,
                   Eigen::VectorXd& p,
                   const Eigen::MatrixXd& X);

#endif
