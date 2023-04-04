//
// Created by raghav on 4/4/23.
//

#ifndef LOADMATRICES_H
#define LOADMATRICES_H

#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

void loadMatrices(const std::vector<std::string>& filepaths,
                  std::vector<Eigen::Matrix4d>& matrices) {
    for (const auto& filepath : filepaths) {
        Eigen::Matrix4d matrix;
        std::ifstream file(filepath);
        if (file.is_open()) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    file >> matrix(i, j);
                }
            }
            matrices.push_back(matrix);
            file.close();
        } else {
            std::cerr << "Unable to open the file: " << filepath << std::endl;
        }
    }
}

#endif