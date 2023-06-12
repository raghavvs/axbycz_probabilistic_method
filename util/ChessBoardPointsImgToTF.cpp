//
// Created by Raghavendra N S on 6/8/23.
//

#include <opencv2/opencv.hpp>
#include <fstream>

int main() {
    // Define the number of corners in your chessboard
    cv::Size patternSize(7, 5);

    std::vector<cv::Point2f> corners;
    std::vector<cv::Point3f> objectPoints;

    // Define real world coordinates for points
    for (int i = 0; i < patternSize.height; i++) {
        for (int j = 0; j < patternSize.width; j++) {
            objectPoints.push_back(cv::Point3f(j, i, 0));
        }
    }

    std::ofstream file("data/CamToBoardTF_new.txt");
    std::ofstream poseFile("data/CamToBoardPose.txt");

    // Loop over all images
    for (int i = 0; i < 20; i++) {
        // Load image
        cv::Mat img = cv::imread("data/nus_eye_on_hand/" + std::to_string(i) + "_img_corner.png");

        // Find chessboard corners
        bool found = cv::findChessboardCorners(img, patternSize, corners);

        if (found) {
            // Compute homography
            cv::Mat H = cv::findHomography(corners, objectPoints);

            // Write homography to file
            file << H << "\n";

            // Normalization to ensure that ||c1|| = 1
            double norm = sqrt(H.at<double>(0,0)*H.at<double>(0,0) +
                               H.at<double>(1,0)*H.at<double>(1,0) +
                               H.at<double>(2,0)*H.at<double>(2,0));
            H /= norm;
            cv::Mat c1  = H.col(0);
            cv::Mat c2  = H.col(1);
            cv::Mat c3 = c1.cross(c2);
            cv::Mat tvec = H.col(2);
            cv::Mat R(3, 3, CV_64F);
            for (int i = 0; i < 3; i++)
            {
                R.at<double>(i,0) = c1.at<double>(i,0);
                R.at<double>(i,1) = c2.at<double>(i,0);
                R.at<double>(i,2) = c3.at<double>(i,0);
            }

            // Construct the 4x4 pose matrix
            cv::Mat pose = cv::Mat::eye(4, 4, CV_64F);
            R.copyTo(pose(cv::Rect(0, 0, 3, 3)));
            tvec.copyTo(pose(cv::Rect(3, 0, 1, 3)));

            // Write pose matrix to file
            poseFile << pose << "\n";
        }
    }

    file.close();
    poseFile.close();

    return 0;
}