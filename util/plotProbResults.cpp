#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <matplotlib> 

namespace plt = matplotlibcpp;

void plotProbResults(std::vector<std::vector<std::vector<double>>> error_1, std::vector<std::vector<std::vector<double>>> error_2, std::vector<double> point, std::string opt) 
{
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

    } else if (opt == "lineplot")
     {
        // Compute the average of the error matrices
        std::vector<std::vector<double>> Err1_Avg(error_1.size(), std::vector<double>(error_1[0][0].size(), 0.0));
        std::vector<std::vector<double>> Err2_Avg(error_2.size(), std::vector<double>(error_2[0][0].size(), 0.0));
        for (size_t i = 0; i < error_1.size(); i++) {
            for (size_t j = 0; j < error_1[0][0].size(); j++) {
                double sum1 = 0.0;
                double sum2 = 0.0;
                for (size_t k = 0; k < error_1[0].size(); k++) {
                    sum1 += error_1[i][k][j];
                    sum2 += error_2[i][k][j];
                }
                Err1_Avg[i][j] = sum1 / error_1[0].size();
                Err2_Avg[i][j] = sum2 / error_2[0].size();
                }
            }
        }
        else if (opt == "lineplot")
        {
        // Compute average values
        VectorXd Err1_Avg = error_1.rowwise().mean();
        VectorXd Err2_Avg = error_2.rowwise().mean();
        
        // Plot line graphs
        plt::subplot(3, 1, 1);
        plt::plot(point, Err1_Avg.segment(0, 3), "o");
        plt::hold(true);
        plt::plot(point, Err2_Avg.segment(0, 3), "*");
        plt::plot(point, Err1_Avg.segment(0, 3), "b");
        plt::plot(point, Err2_Avg.segment(0, 3), "r");
        plt::ylabel("$\\bf E_{R_{X}}$");
        plt::grid(true);
        
        plt::subplot(3, 1, 2);
        plt::plot(point, Err1_Avg.segment(3, 3), "o");
        plt::hold(true);
        plt::plot(point, Err2_Avg.segment(3, 3), "*");
        plt::plot(point, Err1_Avg.segment(3, 3), "b");
        plt::plot(point, Err2_Avg.segment(3, 3), "r");
        plt::ylabel("$\\bf E_{R_{Y}}$");
        plt::grid(true);
        
        plt::subplot(3, 1, 3);
        plt::plot(point, Err1_Avg.segment(6, 3), "o");
        plt::hold(true);
        plt::plot(point, Err2_Avg.segment(6, 3), "*");
        plt::plot(point, Err1_Avg.segment(6, 3), "b");
        plt::plot(point, Err2_Avg.segment(6, 3), "r");
        plt::ylabel("$\\bf E_{R_{Z}}$");
        plt::xlabel("$\\bf Time (sec)$");
        plt::grid(true);
        
        plt::figure();
        plt::subplot(3, 1, 1);
        plt::plot(point, Err1_Avg.segment(9, 3), "o");
        plt::hold(true);
        plt::plot(point, Err2_Avg.segment(9, 3), "*");
        plt::plot(point, Err1_Avg.segment(9, 3), "b");
        plt::plot(point, Err2_Avg.segment(9, 3), "r");
        plt::ylabel("$\\bf E_{t_{X}}$");
        plt::grid(true);
        plt::legend("$Prob1$", "$Prob2$", "fontsize", 12);
        }
}