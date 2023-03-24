// initialize array/vector of matrices - three ways:

#include <iostream>
#include <vector>
#include <array>
#include <eigen3/Eigen/Dense>

int main() {
    /*int len = 10;

    std::vector<Eigen::MatrixXd> C(len);
    C[0] = Eigen::Matrix4d::Random();

    Eigen::MatrixXd D[2];
    D[0] = Eigen::Matrix4d::Random();

    std::array<Eigen::Matrix4d, 2> E;
    D[0] = Eigen::Matrix4d::Random();

    std::cout << "C: \n" << C[0] << std::endl;
    std::cout << "D: \n" << D[0] << std::endl;
    std::cout << "E: \n" << D[0] << std::endl;

    std::cout << "C.size() = " << C.size() << std::endl;
    std::cout << "C.size() = " << C[0].size() << std::endl;*/

    Eigen::MatrixXd m (4,4);
    m << 1, 2, 3, 4,
            5, 6, 7, 8,
            9,10,11,12,
            13,14,15,16;

    std::cout << "Block in the middle" << std::endl;
    std::cout << m.block<2,2> (1,1) << std::endl << std::endl;
    std::cout << m.block<3,1> (0,3) << std::endl << std::endl;
    std::cout << m.block<4,1> (0,3) << std::endl << std::endl;

    return 0;
}

/*
This code uses the .block<p,q>(i,j) function to extract a block of size pxq starting
from row i and column j from the matrix RHS. In this case, it extracts a block of size
3x1 starting from row 0 and column 0. The resulting submatrix is a column vector
containing the first three elements of the first column of RHS.

 using your own words and knowledge, convert the matlab code below to C++:

 using your own words and knowledge, convert the matlab code below to C++,
 and keep in mind to use Eigen:: and std:: explicitly wherever required, also this
 is a continuation of the above code:

 */
