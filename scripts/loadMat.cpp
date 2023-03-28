//
// Created by Raghavendra N S on 3/17/23.
//

// Example code snippet for reading matrices from a text file

#include <iostream>
#include <fstream>
#include <vector>

int main() {
    std::ifstream file("matrix.txt");
    std::vector<std::vector<int>> matrix;
    int rows = 0;
    int cols = 0;
    int num;
    while (file >> num) {
        if (cols == 0) {
            matrix.push_back(std::vector<int>());
            rows++;
        }
        matrix[rows - 1].push_back(num);
        cols++;
        if (cols == 3) {
            cols = 0;
        }
    }
    file.close();
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}