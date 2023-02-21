#include <iostream>
#include <eigen3/Eigen/Dense>
#include <random>
#include <algorithm>
#include <chrono>

using namespace std;
using namespace Eigen;


default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
uniform_real_distribution<double> distribution(0.0, 1.0);

MatrixXd scrambleData(MatrixXd M, double s_rate) {
    int n = M.dimension(2);
    vector<int> M_index(n);
    for (int i = 0; i < n; i++) {
        M_index[i] = i;
    }
    for (int i = 0; i < n; i++) {
        if (distribution(generator) <= 0.01 * s_rate) {
            int index = rand() % n;
            swap(M_index[i], M_index[index]);
        }
    }
    MatrixXd M_perm = MatrixXd::Zero(4, 4, n);
    for (int i = 0; i < n; i++) {
        M_perm.slice(i) = M.slice(M_index[i]);
    }
    return M_perm;
}
