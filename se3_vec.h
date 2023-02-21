#include <eigen3/Eigen/Core>

// Function to vectorize or hat an element of se(3)
Eigen::Matrix<double, 6, 1> se3_vec(Eigen::Matrix4d X) {
    Eigen::Matrix<double, 6, 1> g;
    if (X.cols() == 4) { // If input is skew-sym change to vector
        g << -X(1, 2), X(0, 2), -X(0, 1), X(0, 3), X(1, 3), X(2, 3);
    } else { // If input is vector change to skew-sym
        g << 0, -X(2), X(1), X(3), X(2), -X(0), X(4), -X(1), -X(4), 0, X(5), X(0), -X(3), -X(5), 0, 0;
    }
    return g;
}