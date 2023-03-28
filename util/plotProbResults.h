#ifndef PLOTPROBRESULTS_h
#define PLOTPROBRESULTS_h

#include <Eigen/Dense>
#include <string>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

void plotProbResults(std::vector<std::vector<std::vector<double>>> error_1,
                     std::vector<std::vector<std::vector<double>>> error_2,
                     std::vector<double> point,
                     std::string opt){
    if (opt == "boxplot") {

        // Transpose the error matrices
        std::vector<std::vector<std::vector<double>>> error_1_perm(error_1[0].size(), std::vector<std::vector<double>>(error_1.size(), std::vector<double>(error_1[0][0].size())));
        std::vector<std::vector<std::vector<double>>> error_2_perm(error_2[0].size(), std::vector<std::vector<double>>(error_2.size(), std::vector<double>(error_2[0][0].size())));
        for (size_t i = 0; i < error_1_perm.size(); i++) {
            for (size_t j = 0; j < error_1_perm[0].size(); j++) {
                for (size_t k = 0; k < error_1_perm[0][0].size(); k++) {
                    error_1_perm[i][j][k] = error_1[j][i][k];
                    error_2_perm[i][j][k] = error_2[j][i][k];
                }
            }
        }

        // Plot boxplots for each slice of the permuted error matrices
        for (size_t i = 0; i < error_1_perm[0][0].size(); i++) {
            std::vector<std::vector<double>> data1;
            std::vector<std::vector<double>> data2;
            for (size_t j = 0; j < error_1_perm.size(); j++) {
                data1.push_back(error_1_perm[j][i]);
                data2.push_back(error_2_perm[j][i]);
            }

            plt::figure();
            plt::boxplot(data1, point);
        }
    } else if (opt == "lineplot") {
        // Calculate averages
        Eigen::MatrixXd Err1_Avg = Eigen::MatrixXd::Zero(10, 6);
        Eigen::MatrixXd Err2_Avg = Eigen::MatrixXd::Zero(10, 6);

        for (int i = 0; i < error_1.size(); ++i) {
            for (int j = 0; j < error_1[i].size(); ++j) {
                for (int k = 0; k < error_1[i][j].size(); ++k) {
                    Err1_Avg(i, j) += error_1[i][j][k];
                    Err2_Avg(i, j) += error_2[i][j][k];
                }
                Err1_Avg(i, j) /= error_1[i][j].size();
                Err2_Avg(i, j) /= error_2[i][j].size();
            }
        }

        std::vector<double> point_vec(point.data(), point.data() + point.size());
        std::vector<double> Err1_Avg_row0(Err1_Avg.row(0).data(), Err1_Avg.row(0).data() + Err1_Avg.row(0).size());
        std::vector<double> Err2_Avg_row0(Err2_Avg.row(0).data(), Err2_Avg.row(0).data() + Err2_Avg.row(0).size());
        std::vector<double> Err1_Avg_row1(Err1_Avg.row(1).data(), Err1_Avg.row(1).data() + Err1_Avg.row(1).size());
        std::vector<double> Err2_Avg_row1(Err2_Avg.row(1).data(), Err2_Avg.row(1).data() + Err2_Avg.row(1).size());
        std::vector<double> Err1_Avg_row2(Err1_Avg.row(2).data(), Err1_Avg.row(2).data() + Err1_Avg.row(2).size());
        std::vector<double> Err2_Avg_row2(Err2_Avg.row(2).data(), Err2_Avg.row(2).data() + Err2_Avg.row(2).size());
        std::vector<double> Err1_Avg_row3(Err1_Avg.row(3).data(), Err1_Avg.row(3).data() + Err1_Avg.row(3).size());
        std::vector<double> Err2_Avg_row3(Err2_Avg.row(3).data(), Err2_Avg.row(3).data() + Err2_Avg.row(3).size());
        std::vector<double> Err1_Avg_row4(Err1_Avg.row(4).data(), Err1_Avg.row(4).data() + Err1_Avg.row(4).size());
        std::vector<double> Err2_Avg_row4(Err2_Avg.row(4).data(), Err2_Avg.row(4).data() + Err2_Avg.row(4).size());
        std::vector<double> Err1_Avg_row5(Err1_Avg.row(5).data(), Err1_Avg.row(5).data() + Err1_Avg.row(5).size());
        std::vector<double> Err2_Avg_row5(Err2_Avg.row(5).data(), Err2_Avg.row(5).data() + Err2_Avg.row(5).size());

        plt::figure();
        plt::subplot(3, 1, 1);
        plt::plot(point_vec, Err1_Avg_row0, "o");
        plt::plot(point_vec, Err2_Avg_row0, "*");
        plt::ylabel("$\\bf E_{R_{X}}$", {{"fontsize", "14"}, {"interpreter", "latex"}});
        plt::legend();

        plt::subplot(3, 1, 2);
        plt::plot(point_vec, Err1_Avg_row1, "o");
        plt::plot(point_vec, Err2_Avg_row1, "*");
        plt::ylabel("$\\bf E_{R_{Y}}$", {{"fontsize", "14"}, {"interpreter", "latex"}});
        plt::legend();

        plt::subplot(3, 1, 3);
        plt::plot(point_vec, Err1_Avg_row2, "o");
        plt::plot(point_vec, Err1_Avg_row2, "o");
        plt::plot(point_vec, Err2_Avg_row2, "*");
        plt::ylabel("$\\bf E_{R_{Z}}$", {{"fontsize", "14"}, {"interpreter", "latex"}});
        plt::legend();

        plt::figure();
        plt::subplot(3, 1, 1);
        plt::plot(point_vec, Err1_Avg_row3, "o");
        plt::plot(point_vec, Err2_Avg_row3, "*");
        plt::ylabel("$\\bf E_{t_{Z}}$", {{"fontsize", "14"}, {"interpreter", "latex"}});
        plt::legend();

        // Create plot
        plt::subplot(3, 1, 2);
        plt::plot(point_vec, Err1_Avg_row4, "o");
        plt::plot(point_vec, Err2_Avg_row4, "*");
        plt::plot(point_vec, Err1_Avg_row4, "b");
        plt::plot(point_vec, Err2_Avg_row4, "r");
        plt::ylabel("$\\bf E_{t_{Z}}$", {{"fontsize", "14"}, {"interpreter", "latex"}});
        plt::legend();


        // Create second plot
        plt::subplot(3, 1, 3);
        plt::plot(point_vec, Err1_Avg_row5, "o");
        plt::plot(point_vec, Err2_Avg_row5, "*");
        plt::plot(point_vec, Err1_Avg_row5, "b");
        plt::plot(point_vec, Err2_Avg_row5, "r");
        plt::ylabel("$\\bf E_{t_{Z}}$", {{"fontsize", "14"}, {"interpreter", "latex"}});
        plt::legend();

        // Show plot
        plt::show();
    }
}

#endif