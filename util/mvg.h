/*
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

/*
DESCRIPTION:

The program generates a matrix of multivariate normal samples using the 
Cholesky decomposition method. It takes as input a mean vector, a covariance 
matrix, the number of samples to be generated, and two output matrices to 
hold the generated samples and the lower triangular matrix of the Cholesky 
decomposition of the input covariance matrix. The program first checks that 
the input covariance matrix is positive definite, symmetric and square, and 
that the number of samples requested is a positive integer. It then generates 
the requested number of samples by drawing from a standard normal distribution 
and transforming them using the Cholesky decomposition of the covariance matrix. 
The generated samples and the Cholesky factor are returned as output. 
*/

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#include <random>

void mvg(const Eigen::VectorXd& mu, const Eigen::MatrixXd& Sigma, int N, Eigen::MatrixXd& y, Eigen::MatrixXd& R) {
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
    Eigen::MatrixXd L = Sigma.llt().matrixL();
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
        Eigen::VectorXd r(m);
        for (int j = 0; j < m; j++) {
            r(j) = distribution(generator);
        }
        y.col(i) = L * r + mu;
    }
    R = L;
}

// TEST CASE

int main() {
    Eigen::VectorXd mu(3);
    mu << 1.0, 2.0, 3.0;
    Eigen::MatrixXd Sigma(3, 3);
    Sigma << 2.0, 0.3, 0.5, 0.3, 1.0, -0.4, 0.5, -0.4, 4.0;
    int N = 100;
    Eigen::MatrixXd y;
    Eigen::MatrixXd R;
    mvg(mu, Sigma, N, y, R);
    std::cout << "y = " << std::endl << y << std::endl;
    std::cout << "R = " << std::endl << R << std::endl;
    return 0;
}

