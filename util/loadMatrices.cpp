//
// Created by raghav on 4/4/23.
//

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

int main()
{
    std::vector<Eigen::Matrix4d> A1, B1, C1, A2, B2, C2;

    std::string A1_files = {"data/20230418_abb_charuco_10x14/r1_tf.txt"};
    std::string B1_files = {"data/20230418_abb_charuco_10x14/c2b_tf.txt"};
    std::string C1_files = {"data/20230418_abb_charuco_10x14/r2_tf.txt"};
    std::string A2_files = {"data/20230418_abb_charuco_10x14/r1_tf.txt"};
    std::string B2_files = {"data/20230418_abb_charuco_10x14/c2b_tf.txt"};
    std::string C2_files = {"data/20230418_abb_charuco_10x14/r2_tf.txt"};

    loadMatrices(A1_files, A1);
    loadMatrices(B1_files, B1);
    loadMatrices(C1_files, C1);
    loadMatrices(A2_files, A2);
    loadMatrices(B2_files, B2);
    loadMatrices(C2_files, C2);

    size_t num_matrices = A1.size();
    std::cout << "Number of matrices in A1: " << num_matrices << std::endl;
    std::cout << "A1[0]: " << A1[0] << std::endl;
    std::cout << "A1[1]: " << A1[1] << std::endl;
}
