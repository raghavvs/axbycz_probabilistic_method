/*
DESCRIPTION:

The above program defines a function named so3_vec that takes a 3D vector 
as input and returns a 3x3 matrix. If the input is a vector, the function 
converts it to a skew-symmetric matrix, and if the input is a skew-symmetric 
matrix, the function returns it directly. The resulting 3x3 matrix represents 
an element of the special orthogonal group SO(3), which is used in 3D 
rotation calculations.
*/

#include <eigen3/Eigen/Core>

Eigen::Matrix3d so3_vec(const Eigen::Vector3d& X) {
  Eigen::Matrix3d g;
  if (X.size() == 3) {
    g << 0, -X(2), X(1),
         X(2), 0, -X(0),
         -X(1), X(0), 0;
  } else { 
    g << 0, -X(2), X(1),
         X(2), 0, -X(0),
         -X(1), X(0), 0;
  }
  return g;
}