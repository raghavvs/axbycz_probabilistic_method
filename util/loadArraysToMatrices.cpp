//
// Created by raghav on 11/4/23.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <Eigen/Dense>

void loadArraysToMatrices(const std::vector<std::string>& filepaths,
                          std::vector<Eigen::Matrix4d>& matrices) {
    for (const auto& filepath : filepaths) {
        std::ifstream file(filepath);
        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line)) {
                Eigen::Matrix4d matrix;
                std::istringstream iss(line);
                std::string cell;
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        if (std::getline(iss, cell, ',')) {
                            matrix(i, j) = std::stod(cell);
                        } else {
                            std::cerr << "Error reading matrix element at (" << i << ", " << j << ")" << std::endl;
                            break;
                        }
                    }
                }
                matrices.push_back(matrix);
            }
            file.close();
        } else {
            std::cerr << "Unable to open the file: " << filepath << std::endl;
        }
    }
}

int main() {
    // Replace these file paths with your actual data file paths
    std::vector<std::string> filepaths = {"data/transform_ABC_unified.txt"};
    std::vector<Eigen::Matrix4d> matrices;

    loadArraysToMatrices(filepaths, matrices);

    std::cout << "Loaded " << matrices.size() << " matrices:" << std::endl;
    /*for (const auto& matrix : matrices) {
        std::cout << matrix << std::endl << std::endl;
    }*/

    return 0;
}

