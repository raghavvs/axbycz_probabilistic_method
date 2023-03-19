// initialize array/vector of matrices - three ways:

#include <iostream>
#include <vector>
#include <array>
#include <eigen3/Eigen/Dense>

int main() {
    int len = 10;

    std::vector<Eigen::MatrixXd> C(len);
    C[0] = Eigen::Matrix4d::Random();

    Eigen::MatrixXd D[2];
    D[0] = Eigen::Matrix4d::Random();

    std::array<Eigen::Matrix4d, 2> E;
    D[0] = Eigen::Matrix4d::Random();

    std::cout << "C: \n" << C[0] << std::endl;
    std::cout << "D: \n" << D[0] << std::endl;
    std::cout << "E: \n" << D[0] << std::endl;

    return 0;
}