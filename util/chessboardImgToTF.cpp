//
// Created by Raghavendra N S on 6/8/23.
//

#include <opencv2/opencv.hpp>
#include <iostream>

int main() {
    // Load image
    cv::Mat img = cv::imread("undistorted.png");

    // Define chessboard corners in image coordinates (pixels)
    std::vector<cv::Point2f> img_points;
    img_points.push_back(cv::Point2f(0, 0));
    img_points.push_back(cv::Point2f(img.cols, 0));
    img_points.push_back(cv::Point2f(img.cols, img.rows));
    img_points.push_back(cv::Point2f(0, img.rows));

    // Define chessboard corners in real world coordinates (meters)
    // Change real_points to match your chessboard size (in this case it is 200 X 150 mm)
    std::vector<cv::Point2f> real_points;
    real_points.push_back(cv::Point2f(0, 0));
    real_points.push_back(cv::Point2f(0.2, 0));
    real_points.push_back(cv::Point2f(0.2, 0.15));
    real_points.push_back(cv::Point2f(0, 0.15));

    // Find homography matrix
    cv::Mat H = cv::findHomography(img_points, real_points);

    // Define some data points in image coordinates (pixels)
    std::vector<cv::Point2f> data_points;
    data_points.push_back(cv::Point2f(100, 100));
    data_points.push_back(cv::Point2f(200, 200));
    data_points.push_back(cv::Point2f(300, 300));

    // Apply homography matrix to data points
    std::vector<cv::Point2f> data_points_transformed;
    cv::perspectiveTransform(data_points, data_points_transformed, H);

    // Print data points before and after transformation
    std::cout << "Data points before transformation:\n";
    for (const auto& point : data_points)
        std::cout << point << "\n";

    std::cout << "Data points after transformation:\n";
    for (const auto& point : data_points_transformed)
        std::cout << point << "\n";

    // Draw data points on image
    for (const auto& point : data_points)
        cv::circle(img, point, 5, cv::Scalar(0, 255, 0), -1);

    // Save image with data points
    cv::imwrite("undistorted_withPoints.png", img);

    // Apply homography matrix to image
    cv::Mat new_img;
    cv::warpPerspective(img, new_img, H, cv::Size(int(0.3 * img.cols), int(0.3 * img.rows)));

    // Save transformed image
    cv::imwrite("new_undistorted.png", new_img);

    // Print homography matrix
    std::cout << "Homography matrix:\n" << H << "\n";

    return 0;
}