//
// Created by Raghavendra N S on 3/17/23.
//

// Example code snippet for reading matrices from a text file

#include <iostream>
#include <fstream>
#include <Eigen/Dense>

int main() {
    std::string file_path = "data/000.txt";
    std::ifstream file(file_path);
    Eigen::Matrix4d matrix;

    if (file.is_open()) {
        std::string line;
        int row = 0;

        while (std::getline(file, line) && row < 4) {
            std::istringstream iss(line);
            double value;
            int col = 0;

            while (iss >> value && col < 4) {
                matrix(row, col) = value;
                col++;
            }

            row++;
        }
    } else {
        std::cerr << "Unable to open the file: " << file_path << std::endl;
        return 1;
    }

    // Print the 4x4 matrix
    std::cout << matrix << std::endl;

    return 0;
}
