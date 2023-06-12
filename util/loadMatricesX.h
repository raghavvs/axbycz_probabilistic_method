//
// Created by Raghavendra N S on 6/12/23.
// For NUS file formatting
//

#ifndef LOADMATRICESX_H
#define LOADMATRICESX_H

#include <Eigen/Dense>
#include <fstream>
#include <vector>

void loadMatricesX(std::vector<Eigen::Matrix4d> &A,
                   std::vector<Eigen::Matrix4d> &B) {
    for (int i = 0; i <= 20; i++) {
        std::string markerFileName = "data/221221_panda_ps_EH/" + std::to_string(i) + "_markerpose.txt";
        std::string robotFileName = "data/221221_panda_ps_EH/" + std::to_string(i) + "_robotpose.txt";

        Eigen::Matrix4d markerMatrix;
        Eigen::Matrix4d robotMatrix;

        std::ifstream markerFile(markerFileName);
        if (markerFile.is_open()) {
            for (int row = 0; row < 4; row++) {
                for (int col = 0; col < 4; col++) {
                    markerFile >> markerMatrix(row, col);
                }
            }
            B.push_back(markerMatrix);
        }
        markerFile.close();

        std::ifstream robotFile(robotFileName);
        if (robotFile.is_open()) {
            for (int row = 0; row < 4; row++) {
                for (int col = 0; col < 4; col++) {
                    robotFile >> robotMatrix(row, col);
                }
            }
            A.push_back(robotMatrix);
        }
        robotFile.close();
    }
}

#endif