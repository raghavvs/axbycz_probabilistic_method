#include <Eigen/Dense>
#include <string>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

void plotProbResults(const Eigen::MatrixXd &error_1, const Eigen::MatrixXd &error_2, const Eigen::VectorXd &point, const std::string &opt) {
    if (opt == "boxplot") {
        Eigen::MatrixXd error_1_perm = error_1.transpose();
        Eigen::MatrixXd error_2_perm = error_2.transpose();

        for (int i = 0; i < error_1_perm.cols(); ++i) {
            // figure
            // boxplot(error_1_perm.col(i), point);
        }

        for (int i = 0; i < error_2_perm.cols(); ++i) {
            // figure
        }
    } else if (opt == "lineplot") {
        Eigen::RowVectorXd Err1_Avg = error_1.colwise().sum() / error_1.cols();
        Eigen::RowVectorXd Err2_Avg = error_2.colwise().sum() / error_2.cols();

        plt::figure();
        plt::subplot(3, 1, 1);
        plt::plot(point, Err1_Avg.row(0), "o");
        plt::plot(point, Err2_Avg.row(0), "*");
        plt::plot(point, Err1_Avg.row(0), "b");
        plt::plot(point, Err2_Avg.row(0), "r");
        plt::ylabel("$\\bf E_{R_{X}}$", {{"fontsize", 14}, {"interpreter", "latex"}});
        plt::legend({"$Prob1$", "$Prob2$"}, {{"fontsize", 14}, {"interpreter", "latex"}});

        plt::subplot(3, 1, 2);
        plt::plot(point, Err1_Avg.row(1), "o");
        plt::plot(point, Err2_Avg.row(1), "*");
        plt::plot(point, Err1_Avg.row(1), "b");
        plt::plot(point, Err2_Avg.row(1), "r");
        plt::ylabel("$\\bf E_{R_{Y}}$", {{"fontsize", 14}, {"interpreter", "latex"}});
    }
}