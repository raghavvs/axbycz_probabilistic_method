#include <iostream>
#include <Eigen/Dense>
#include <generateABC.h>

using namespace Eigen;

void generateSetsOfABC(int Num, int optPDF, Vector3d Mean, Matrix3d Cov, VectorXd XActual, VectorXd YActual, VectorXd ZActual, 
                       VectorXd& A1, VectorXd& B1, VectorXd& C1, VectorXd& A2, VectorXd& B2, VectorXd& C2, VectorXd& A3, VectorXd& B3, VectorXd& C3) {

    // Generate constant A1, free B1 and C1
    int opt = 1;
    generateABC(Num, opt, optPDF, Mean, Cov, XActual, YActual, ZActual, A1, B1, C1);

    // Generate constant C2, free A2 and B2
    opt = 3;
    generateABC(Num, opt, optPDF, Mean, Cov, XActual, YActual, ZActual, A2, B2, C2);

    // Generate constant B3, free A3 and C3
    opt = 2;
    generateABC(Num, opt, optPDF, Mean, Cov, XActual, YActual, ZActual, A3, B3, C3);
}
