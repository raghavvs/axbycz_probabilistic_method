//
// Created by Raghavendra N S on 3/17/23.
//

// Example code snippet for reading matrices from a text file

#include <fstream>
#include <vector>
#include <Eigen/Dense>

int main() {
    std::ifstream file("data.txt");
    std::vector<Eigen::Matrix4d> MatA;
    for (int k = 0; k < 10; ++k) {
        Eigen::Matrix4d mat;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                file >> mat(i,j);
        MatA.push_back(mat);
    }
}