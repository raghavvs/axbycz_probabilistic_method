#include <eigen3/Eigen/Core>

Eigen::Matrix3d so3Vec(const Eigen::Vector3d& X) {
  Eigen::Matrix3d g;
  if (X.size() == 3) { // If input is skew-sym change to vector
    g << 0, -X(2), X(1),
         X(2), 0, -X(0),
         -X(1), X(0), 0;
  } else { // If input is vector change to skew-sym
    g << 0, -X(2), X(1),
         X(2), 0, -X(0),
         -X(1), X(0), 0;
  }
  return g;
}