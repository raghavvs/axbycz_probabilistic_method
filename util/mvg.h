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

MVG - Multivariate Gaussian random number generator.

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

The function takes in four arguments:

    mu: a vector of means of length m, where m is the dimension of the 
        multivariate normal distribution.
    Sigma: a square matrix of size m x m representing the covariance matrix.
    N: an integer representing the number of samples to generate.
    y: a matrix of size m x N where the generated samples will be stored.

The function first checks that the input arguments meet certain requirements:

    mu must have the same length as the number of rows in Sigma.
    Sigma must be square.
    Sigma must be symmetric.
    Sigma must be positive definite.
    N must be a positive intege
*/

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#include <random>

std::pair<Eigen::VectorXd, Eigen::MatrixXd> mvg(const Eigen::VectorXd& mu, const Eigen::MatrixXd& Sigma, int N) {
    if (mu.size() != Sigma.rows()) {
        std::cerr << "Length(mu) must equal size(Sigma,1)." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (Sigma.rows() != Sigma.cols()) {
        std::cerr << "Sigma must be square." << std::endl;
        exit(EXIT_FAILURE);
    }

    if ((Sigma - Sigma.transpose()).norm() > 1e-15) {
        std::cerr << "Sigma must be symmetric." << std::endl;
        exit(EXIT_FAILURE);
    }

    bool is_spd;
    Eigen::MatrixXd L = Sigma.llt().matrixL();
    is_spd = ((L * L.transpose() - Sigma).norm() < 1e-8);

    if (!is_spd) {
        std::cerr << "Sigma must be positive definite." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (N < 1) {
        std::cerr << "A positive integer number of samples must be requested." << std::endl;
        exit(EXIT_FAILURE);
    }

    int m = mu.size();
    Eigen::MatrixXd y(m, N);
    Eigen::MatrixXd R = L;

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < N; i++) {
        Eigen::VectorXd r(m);
        for (int j = 0; j < m; j++) {
            r(j) = distribution(generator);
        }
        y.col(i) = L * r + mu;
    }

    return {y, R};
}

// TEST CASE

/* int main() {
     // Define mean vector
    Eigen::VectorXd mu(4);
    mu << 1, 2, 3, 4;

    // Define covariance matrix
    Eigen::MatrixXd Sigma(4, 4);
    Sigma << 2, 0.5, -1, 0,
             0.5, 1, 0, -1,
             -1, 0, 3, 0.5,
             0, -1, 0.5, 2;

    int N = 10; // number of samples to be generated

    Eigen::VectorXd y;
    Eigen::MatrixXd R;

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> result = mvg(mu, Sigma, N);
    y = result.first;
    R = result.second;

    std::cout << "Generated samples:\n" << y << std::endl;
    std::cout << "Cholesky factorization of covariance matrix:\n" << R << std::endl;

    return 0;
} */
