//
// Created by Raghavendra N S on 6/8/23.
//

#include <opencv2/opencv.hpp>
#include <fstream>

int main() {
    // Define the number of corners in your chessboard
    cv::Size patternSize(7, 5);  // for example, a 9x6 chessboard

    std::vector<cv::Point2f> corners;
    std::vector<cv::Point3f> objectPoints;

    // Define real world coordinates for points
    for (int i = 0; i < patternSize.height; i++) {
        for (int j = 0; j < patternSize.width; j++) {
            objectPoints.push_back(cv::Point3f(j, i, 0));
        }
    }

    std::ofstream file("CamToBoardTF.txt");

    // Loop over all images
    for (int i = 0; i < 10; i++) {
        // Load image
        cv::Mat img = cv::imread("data/nus_eye_on_hand/" + std::to_string(i) + "+_img_corner.png");

        // Find chessboard corners
        bool found = cv::findChessboardCorners(img, patternSize, corners);

        if (found) {
            // Compute homography
            cv::Mat H = cv::findHomography(corners, objectPoints);

            // Write homography to file
            file << "Homography for image " << i << ":\n";
            file << H << "\n";
        }
    }

    file.close();

    return 0;
}