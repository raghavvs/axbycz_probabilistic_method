//
// Created by raghav on 4/4/23.
//

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

int main()
{
    std::vector<Eigen::Matrix4d> A1, C1, A2, C2;

    std::vector<std::string> A1_files = {"data/000.txt"};
    std::vector<std::string> C1_files = {"data/000_m.txt"};
    std::vector<std::string> A2_files = {"data/001.txt"};
    std::vector<std::string> C2_files = {"data/001_m.txt"};

    loadMatrices(A1_files, A1);
    loadMatrices(C1_files, C1);
    loadMatrices(A2_files, A2);
    loadMatrices(C2_files, C2);

    std::cout << "A1: " << A1[0] << std::endl;
    std::cout << "C1: " << C1[0] << std::endl;
    std::cout << "A2: " << A2[0] << std::endl;
    std::cout << "C2: " << C2[0] << std::endl;
}
