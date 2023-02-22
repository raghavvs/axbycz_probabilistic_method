/*
MVG    Multivariate Gaussian random number generator.

   y = mvg(mu,Sigma,N), where mu is mx1 and Sigma is mxm and SPD, produces 
   an mxN matrix y whose columns are samples from the multivariate 
   Gaussian distribution parameterized by mean mu and covariance Sigma.

   [y,R] = mvg(mu,Sigma,N) also returns the Cholesky factor of the
   covariance matrix Sigma such that Sigma = R'*R.

   See also RAND, RANDN, SPRANDN, SPRANDN, RANDPERM.

   Chad Lieberman, MIT 2008.
   Questions/Comments:  celieber@mit.edu
   $Revision: 1.0.0 $  $Date: 2008/09/01 $
   $Revision: 1.0.1 $  $Date: 2008/09/03 $

   References:  
   [1] I.T. Hernadvolgyi (1998) "Generating random vectors from the 
       multivariate normal distribution."
   Available on-line at http://www.csi.uottawa.ca/~istvan/work.html.

   Acknowledgements:  
       I would like to acknowledge John D'Errico for his helpful comments 
       and suggestions.
*/

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#include <random>

using namespace Eigen;

void mvg(const VectorXd& mu, const MatrixXd& Sigma, int N, MatrixXd& y, MatrixXd& R) {
    if (mu.size() != Sigma.rows()) {
        std::cerr << "Length(mu) must equal size(Sigma,1)." << std::endl;
        return;
    }
    if (Sigma.rows() != Sigma.cols()) {
        std::cerr << "Sigma must be square." << std::endl;
        return;
    }
    if ((Sigma - Sigma.transpose()).norm() > 1e-15) {
        std::cerr << "Sigma must be symmetric." << std::endl;
        return;
    }
    bool is_spd;
    MatrixXd L = Sigma.llt().matrixL();
    is_spd = ((L * L.transpose() - Sigma).norm() < 1e-8);
    if (!is_spd) {
        std::cerr << "Sigma must be positive definite." << std::endl;
        return;
    }
    if (N < 1) {
        std::cerr << "A positive integer number of samples must be requested." << std::endl;
        return;
    }
    int m = mu.size();
    y.resize(m, N);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < N; i++) {
        VectorXd r(m);
        for (int j = 0; j < m; j++) {
            r(j) = distribution(generator);
        }
        y.col(i) = L * r + mu;
    }
    R = L;
}

int main() {
    VectorXd mu(3);
    mu << 1.0, 2.0, 3.0;
    MatrixXd Sigma(3, 3);
    Sigma << 2.0, 0.3, 0.5, 0.3, 1.0, -0.4, 0.5, -0.4, 4.0;
    int N = 100;
    MatrixXd y;
    MatrixXd R;
    mvg(mu, Sigma, N, y, R);
    std::cout << "y = " << std::endl << y << std::endl;
    std::cout << "R = " << std::endl << R << std::endl;
    return 0;
}
