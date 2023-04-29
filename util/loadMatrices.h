//
// Created by raghav on 4/4/23.
//

#ifndef LOADMATRICES_H
#define LOADMATRICES_H

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void loadMatrices(const std::string& filepath,
                  std::vector<Eigen::Matrix4d>& matrices) {
    Eigen::Matrix4d matrix;
    std::ifstream file(filepath);
    if (file.is_open()) {
        while (!file.eof()) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    file >> matrix(i, j);
                }
            }
            matrices.push_back(matrix);
        }
        file.close();
    } else {
        std::cerr << "Unable to open the file: " << filepath << std::endl;
    }
}

#endif