//
// Created by raghav on 4/4/23.
//
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <vector>
#include <filesystem>

std::string find_file(const std::string &filename) {
    for (const auto &entry : std::filesystem::recursive_directory_iterator(std::filesystem::current_path())) {
        if (entry.path().filename() == filename) {
            return entry.path();
        }
    }
    return "";
}

int main() {
    std::string file_path = find_file("tf_echo.txt");

    if (file_path.empty()) {
        std::cerr << "Error: unable to find file 'tf_echo.txt'" << std::endl;
        return 1;
    }

    std::ifstream file(file_path);
    // ...
}
