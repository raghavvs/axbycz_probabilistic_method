//
// Created by raghav on 4/4/23.
//
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <vector>

int main() {
    std::ifstream file("tf_echo.txt");

    if (!file.is_open()) {
        std::cerr << "Error: unable to open file 'tf_echo.txt'" << std::endl;
        return 1;
    }

    std::string data((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    std::regex translation_pattern("Translation: \\[(.+?)\\]");
    std::regex rotation_pattern("Rotation: in Quaternion \\[(.+?)\\]");

    std::vector<std::string> translations;
    std::vector<std::string> rotations;

    for (std::sregex_iterator i = std::sregex_iterator(data.begin(), data.end(), translation_pattern); i != std::sregex_iterator(); ++i) {
        translations.push_back((*i)[1].str());
    }

    for (std::sregex_iterator i = std::sregex_iterator(data.begin(), data.end(), rotation_pattern); i != std::sregex_iterator(); ++i) {
        rotations.push_back((*i)[1].str());
    }

    for (size_t i = 0; i < translations.size(); ++i) {
        std::string translation = translations[i];
        std::string rotation = rotations[i];

        std::cout << "Translation: [";
        size_t pos = 0;
        while ((pos = translation.find(',')) != std::string::npos) {
            std::cout << stof(translation.substr(0, pos)) << ", ";
            translation.erase(0, pos + 1);
        }
        std::cout << stof(translation) << "]" << std::endl;

        std::cout << "Rotation: [";
        pos = 0;
        while ((pos = rotation.find(',')) != std::string::npos) {
            std::cout << stof(rotation.substr(0, pos)) << ", ";
            rotation.erase(0, pos + 1);
        }
        std::cout << stof(rotation) << "]" << std::endl;
    }

    return 0;
}